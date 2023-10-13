import datetime
import os
import appdirs
import json
import re
import logging
import Bio.Entrez
import lxml.etree
import colorama
import pyparsing
import urllib
import sys


PROJECTNAME = 'paperboy'
CONFFILE = '{}.conf'.format(PROJECTNAME)
EMAIL = 'sten.heinze@gmail.de'
DATEFORMAT_INTERNAL = '%Y-%m-%d'
DATEFORMAT_ENTREZ = '%Y/%m/%d'
DATEFORMAT_USER = '%Y-%m-%d'

# run on import

# use min width so that output is not shortened so much to become unreadable
MIN_TERMINAL_WIDTH = 110
try:
    MAX_LEN = max(MIN_TERMINAL_WIDTH, os.get_terminal_size()[0])
except OSError:
    MAX_LEN = sys.maxsize

# https://stackoverflow.com/questions/2186919/getting-correct-string-length-in-python-for-strings-with-ansi-color-codes
ESC = pyparsing.Literal('\x1b')
integer = pyparsing.Word(pyparsing.nums)
escapeSeq = pyparsing.Combine(
    ESC + '[' + pyparsing.Optional(pyparsing.delimitedList(integer, ';')) + pyparsing.oneOf(list(pyparsing.alphas)))


def nonAnsiString(s): return pyparsing.Suppress(escapeSeq).transformString(s)


colorama.init()


# https://gist.github.com/setaou/ff98e82a9ce68f4c2b8637406b4620d1
class JSONDecoder2(json.JSONDecoder):
    date_regex = re.compile(r'(\d{1,4}[-/]\d{1,2}[-/]\d{1,2})')

    def __init__(self, *args, **kwargs):
        logging.debug('JSONDecoder2.__init__')
        json.JSONDecoder.__init__(
            self, object_hook=self.object_hook_func, *args, **kwargs)
        self.parse_string = self.new_scanstring
        # Use python version, the C version does not use the new parse_string
        self.scan_once = json.scanner.py_make_scanner(self)

    @classmethod
    def new_scanstring(cls, s, end, strict=True):
        (s, end) = json.decoder.scanstring(s, end, strict)
        logging.debug('JSONDecoder2.new_scanstring: {}'.format(s))
        if cls.date_regex.match(s):
            logging.debug('JSONDecoder2.new_scanstring: date_regex matches')
            return (datetime.date.fromisoformat(s), end)
        else:
            return (s, end)

    def object_hook_func(self, obj):
        if '__type__' in obj and obj['__type__'] == 'Journal':
            return Journal(obj['__data__'])
        return obj


def custom_encoder(obj):
    if isinstance(obj, datetime.date):
        return obj.isoformat()
    if isinstance(obj, Journal):
        return {'__type__': 'Journal', '__data__': obj.medline_format()}
    raise TypeError('Type {} not JSON serializable'.format(type(obj)))


class Config:
    def __init__(self):
        self.last_check_date = datetime.date.today() - datetime.timedelta(days=4)
        self.email = EMAIL
        self.pubmed_journals_url = "https://ftp.ncbi.nih.gov/pubmed/J_Medline.txt"
        self.pubmed_journals = []
        self.pubmed_journals_last_update = datetime.date.min
        self.pubmed_journal_update_interval_days = 30
        self.journals = []
        # dict {journal_string : [article_id, ..], ..}
        self.journals_last_article_id_lists = dict()

        self.load_from_file()

    def get_member_dict(self):
        member_dict = {attr: getattr(self, attr) for attr in dir(self)
                       if not callable(getattr(self, attr)) and not attr.startswith("__")}
        logging.debug('Config.get_member_dict: {}'.format(member_dict))
        return member_dict

    def set_members(self, member_dict):
        members = [attr for attr in dir(self)
                   if not callable(getattr(self, attr)) and not attr.startswith("__")]
        for m in members:
            if m in member_dict:
                setattr(self, m, member_dict[m])
        logging.debug('Config.set_members: {}'.format(member_dict))

    def load_from_file(self):
        cfg_file = os.path.join(appdirs.user_config_dir(PROJECTNAME), CONFFILE)

        try:
            with open(cfg_file, 'r') as f:
                self.set_members(json.load(f, cls=JSONDecoder2))
                logging.debug(
                    'Config.load_from_file: loading from {}'.format(cfg_file))
        except OSError:
            logging.debug(
                'Config.load_from_file: could not read {}, using defaults'.format(cfg_file))

    def save_to_file(self):
        folder = appdirs.user_config_dir(PROJECTNAME)
        os.makedirs(folder, exist_ok=True)

        cfg_file = os.path.join(folder, CONFFILE)
        logging.debug('Config.save_to_file: saving to {}'.format(cfg_file))

        with open(cfg_file, 'w') as f:
            json.dump(self.get_member_dict(), f,
                      default=custom_encoder, indent=4, sort_keys=True)


class Article:
    def __init__(self, article_xml_tree):
        xpath_doi = './PubmedData/ArticleIdList/ArticleId[@IdType="doi"]'
        xpath_pmid = './PubmedData/ArticleIdList/ArticleId[@IdType="pubmed"]'
        xpath_title = './MedlineCitation/Article/ArticleTitle'
        xpath_pub_type = './MedlineCitation/Article/PublicationTypeList/PublicationType'
        xpath_journal_abbrev = './MedlineCitation/Article/Journal/ISOAbbreviation'
        xpath_journal_vol = './MedlineCitation/Article/Journal/JournalIssue/Volume'
        xpath_journal_issue = './MedlineCitation/Article/Journal/JournalIssue/Issue'
        xpath_journal_y = './MedlineCitation/Article/Journal/JournalIssue/PubDate/Year'
        xpath_journal_m = './MedlineCitation/Article/Journal/JournalIssue/PubDate/Month'
        xpath_journal_d = './MedlineCitation/Article/Journal/JournalIssue/PubDate/Day'
        xpath_authorlist = './MedlineCitation/Article/AuthorList'
        xpath_author_lastname = './Author/LastName'
        # xpath_author_firstname = './Author/ForeName'
        # xpath_author_middle = './Author/Initials'
        xpath_author_collective_name = './Author/CollectiveName'

        self.missing_data = []
        self.doi = 'https://doi.org/{}'.format(
            self.extract_first_data(article_xml_tree, 'DOI', xpath_doi))
        self.pmid = self.extract_first_data(
            article_xml_tree, 'PMID', xpath_pmid)
        self.title = self.extract_first_data(
            article_xml_tree, 'Title', xpath_title)
        self.pub_type_list = self.extract_all_data_list(
            article_xml_tree, 'PubType', xpath_pub_type)
        self.journal_abbrev = self.extract_first_data(
            article_xml_tree, 'JournalAbbrev', xpath_journal_abbrev)
        self.journal_vol = self.extract_first_data(
            article_xml_tree, 'JournalVol', xpath_journal_vol)
        self.journal_issue = self.extract_first_data(
            article_xml_tree, 'JournalIssue', xpath_journal_issue)

        y = self.extract_first_data(
            article_xml_tree, 'JournalYear', xpath_journal_y)
        m = self.extract_first_data(
            article_xml_tree, 'JournalMonth', xpath_journal_m)
        d = self.extract_first_data(
            article_xml_tree, 'JournalDay', xpath_journal_d)
        logging.debug('Article.__init__: y={} m={} d={}'.format(y, m, d))
        # some journals do not use day or month field,
        # set to 1 to avoid breaking the date parsing
        if m is None:
            m = 1
        if d is None:
            d = 1
        try:
            self.journal_date = datetime.datetime.strptime(
                '{}-{}-{}'.format(y, m, d), '%Y-%b-%d').date()
        except ValueError:
            try:
                self.journal_date = datetime.strptime(
                    '{}-{}-{}'.format(y, m, d), '%Y-%m-%d').date()
            except ValueError:
                self.missing_data.append('JournalDate')

        authorlist_xml_tree_list = article_xml_tree.xpath(xpath_authorlist)
        if not len(authorlist_xml_tree_list) == 1:
            self.missing_data.append('AuthorList')
            self.author_lastname_list = []
        else:
            authorlist_xml_tree = authorlist_xml_tree_list[0]
            self.author_lastname_list = self.extract_all_data_list(
                authorlist_xml_tree, 'Lastname', xpath_author_lastname)
            if self.author_lastname_list is None:
                self.author_lastname_list = self.extract_all_data_list(
                    authorlist_xml_tree, 'CollectiveName', xpath_author_collective_name)
        logging.debug('Article.__init__: authors={}'.format(
            self.author_lastname_list))

        logging.debug('Article.__init__: {}'.format(self.__repr__()))
        if len(self.missing_data) > 0:
            logging.debug('Missing data for PMID {}: {}'.format(
                self.pmid, ', '.join(self.missing_data)))

    # object representation
    def __repr__(self):
        return 'Article(DOI={}, PMID={}, \'{}\' {}. {}, {}, {}, {}, {})' \
            .format(self.doi, self.pmid, self.title,
                    ', '.join(self.author_lastname_list),
                    self.journal_abbrev, self.journal_vol, self.journal_issue,
                    self.journal_date.strftime(DATEFORMAT_INTERNAL),
                    ', '.join(self.pub_type_list))

    # user readable object information that fits on one terminal line (MAX_LEN)
    def __str__(self):
        # some article have no authors, e.g. news, but are still categorized
        # as 'Journal Article'
        author_lastnames = ''
        if len(self.author_lastname_list) > 0:
            author_lastnames = '{}. '.format(
                ', '.join(self.author_lastname_list))

        # article can have many publication types, remove less informative ones
        pub_types_to_remove = ['Research Support, Non-U.S. Gov\'t',
                               'Research Support, U.S. Gov\'t, Non-P.H.S.',
                               'Research Support, U.S. Gov\'t, P.H.S.',
                               'Research Support, N.I.H., Extramural',
                               'Research Support, N.I.H., Intramural']
        pub_type_list_to_show = [p for p in self.pub_type_list
                                 if p not in pub_types_to_remove]

        journal_doi_pmid = '{} ('.format(self.journal_abbrev) \
            + colorama.Fore.YELLOW + '{}'.format(self.doi) + colorama.Style.RESET_ALL \
            + ') (PMID ' \
            + colorama.Fore.YELLOW + '{}'.format(self.pmid) + colorama.Style.RESET_ALL \
            + ') ({})'.format(', '.join(pub_type_list_to_show))

        title = colorama.Style.BRIGHT + '{}'.format(self.title) + colorama.Style.RESET_ALL + ' '

        result = title + author_lastnames + journal_doi_pmid

        if len(nonAnsiString(result)) > MAX_LEN and len(self.author_lastname_list) > 2:
            logging.debug('Article.__str__: len_str={} > MAX_LEN={}, len_auth={} > 2'.format(
                len(nonAnsiString(result)), MAX_LEN, len(self.author_lastname_list)))
            author_lastnames = ', '.join(self.author_lastname_list[0:1]) \
                + ', ({} more). '.format(len(self.author_lastname_list) - 1)
            result = title + author_lastnames + journal_doi_pmid

        if len(nonAnsiString(result)) > MAX_LEN:
            logging.debug('Article.__str__: len_str={} > MAX_LEN={}'.format(
                len(nonAnsiString(result)), MAX_LEN))
            title_len = len(self.title) + MAX_LEN - \
                len(nonAnsiString(result)) - 2
            title = colorama.Style.BRIGHT \
                + '{}..'.format(self.title[0:title_len].strip()) \
                + colorama.Style.RESET_ALL + ' '
            result = title + author_lastnames + journal_doi_pmid

        return result

    def extract_first_data(self, xml_tree, description, xpath_spec):
        data_list = xml_tree.xpath(xpath_spec)
        if len(data_list) == 1:
            logging.debug('Article.extract_first_data: {}'.format(
                lxml.etree.tostring(data_list[0])))

            # https://stackoverflow.com/questions/4624062/get-all-text-inside-a-tag-in-lxml
            data_string = ''.join(data_list[0].itertext())
            logging.debug('Article.extract_first_data: {}'.format(data_string))
            return data_string
        else:
            self.missing_data.append(description)
            return None

    def extract_all_data_list(self, xml_tree, description, xpath_spec):
        data_list = xml_tree.xpath(xpath_spec)
        if len(data_list) > 0:
            logging.debug('Article.extract_all_data_list: {}'.format(
                [lxml.etree.tostring(d) for d in data_list]))
            data_string_list = [''.join(d.itertext()) for d in data_list]
            logging.debug('Article.extract_all_data_list: {}'.format(
                data_string_list))
            return [d.replace('\r', '').replace('\n', '') for d in data_string_list]
        else:
            self.missing_data.append(description)
            return None


class Journal:
    j_title_key = 'JournalTitle'
    j_medabbr_key = 'MedAbbr'
    j_issn_print_key = 'ISSN (Print)'
    j_issn_online_key = 'ISSN (Online)'
    j_nlmid_key = 'NlmId'

    def __init__(self, medline_journal_data):
        journal_data_split_char = ':'
        self.journal_data_dict = dict()
        for line in medline_journal_data:
            line_parts = line.split(journal_data_split_char, 1)
            if len(line_parts) < 2:
                logging.debug('Journal.__init__: Ignoring line \'{}\''.format(
                    line))
                continue

            data = line_parts[1].strip()
            if data:
                self.journal_data_dict[line_parts[0].strip()] = data

    def medline_format(self):
        return ['{}: {}'.format(k, v) for k, v in self.journal_data_dict.items()]

    def __repr__(self):
        return 'Journal({})'.format(self.journal_data_dict)

    def __str__(self):
        title = self.journal_data_dict[Journal.j_title_key] if Journal.j_title_key in self.journal_data_dict else "<No Title>"
        abbr_part = ' ({})'.format(
            self.journal_data_dict[Journal.j_medabbr_key]) \
            if Journal.j_medabbr_key in self.journal_data_dict else ""
        nlmid_part = ' ({}: '.format(Journal.j_nlmid_key) \
            + colorama.Fore.YELLOW \
            + '{}'.format(self.journal_data_dict[Journal.j_nlmid_key]) \
            + colorama.Style.RESET_ALL + ')' \
            if Journal.j_nlmid_key in self.journal_data_dict else ""

        issn_list = []
        if Journal.j_issn_print_key in self.journal_data_dict:
            issn_list.append(self.journal_data_dict[Journal.j_issn_print_key])
        if Journal.j_issn_online_key in self.journal_data_dict:
            issn_list.append(self.journal_data_dict[Journal.j_issn_online_key])
        issn_part = ""
        if len(issn_list) > 0:
            issn_part = ' (ISSN: ' + ', '.join(issn_list) + ')'

        result = 'Journal: ' + colorama.Style.BRIGHT + \
            '{}'.format(title) + colorama.Style.RESET_ALL + \
            abbr_part + nlmid_part + issn_part

        if len(nonAnsiString(result)) > MAX_LEN:
            logging.debug('Journal.__str__: len_str={} > MAX_LEN={}'.format(
                len(nonAnsiString(result)), MAX_LEN))
            title_len = len(title) + MAX_LEN - len(nonAnsiString(result)) - 2
            result = 'Journal: ' + Style.BRIGHT + \
                '{}..'.format(title[0:title_len]) + Style.RESET_ALL + \
                abbr_part + nlmid_part + issn_part

        return result


class Paperboy:
    # constructor
    def __init__(self):
        self.cfg = Config()
        self.articles = []
        Bio.Entrez.email = self.cfg.email

    def entrez_esearch(self, *args, **kwargs):
        logging.debug('Paperboy.entrez_esearch: {} {}'.format(args, kwargs))

        try:
            handle = Bio.Entrez.esearch(*args, **kwargs)
            record = Bio.Entrez.read(handle)
            handle.close()
            return record
        except (IOError, OSError):
            logging.error(
                'Paperboy.entrez_esearch: I/O Error searching data from NCBI Entrez.')
            sys.exit()

    def entrez_efetch(self, *args, **kwargs):
        logging.debug('Paperboy.entrez_efetch: {} {}'.format(args, kwargs))

        try:
            handle = Bio.Entrez.efetch(*args, **kwargs)
            record = handle.read()
            handle.close()
            return record
        except (IOError, OSError):
            logging.error(
                'Paperboy.entrez_efetch: I/O Error fetching data from NCBI Entrez.')
            sys.exit()

    # load and return all (up to 100000) article IDs for a given journal_nlmid
    # and date range (pubmed does not support time ranges, only date ranges)
    # can include previously found articles from the same date range
    def load_article_ids_from_journal(self, journal_nlmid):
        # check until today (pubmed seems to use inclusive ranges)
        term = journal_nlmid + '[jid] AND ' \
            + self.cfg.last_check_date.strftime(DATEFORMAT_ENTREZ) + ":" \
            + datetime.date.today().strftime(DATEFORMAT_ENTREZ) + '[dp]'
        record = self.entrez_esearch(db='pubmed', term=term)
        logging.debug('Paperboy.load_article_ids_from_journal: {}'.format(
            record))

        article_ids = record['IdList']

        # limit to 100000, otherwise needs retstart loop
        # https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
        count = min(int(record['Count']), 100000)
        retsize = len(article_ids)

        if count > retsize:
            logging.debug('Paperboy.load_article_ids_from_journal: count={} > retsize={}'.format(
                count, retsize))
            record = self.entrez_esearch(
                db='pubmed', term=term, retmax=count, retstart=retsize)
            logging.debug('Paperboy.load_article_ids_from_journal: {}'.format(
                record))

            article_ids.extend(record['IdList'])

        return article_ids

    # load article IDs for all journals
    def update_article_ids(self):
        active_journals = []
        for pubmed_j in self.cfg.pubmed_journals:
            for active_j in self.cfg.journals:
                if active_j == pubmed_j.journal_data_dict[Journal.j_nlmid_key]:
                    active_journals.append(pubmed_j)
        logging.debug('Paperboy.update_article_ids: updating article ids for {} journals'.format(
            len(active_journals)))

        for j in active_journals:
            j_nlmid = j.journal_data_dict[Journal.j_nlmid_key]
            j_abbr = j.journal_data_dict[Journal.j_medabbr_key]
            # load article IDs
            article_ids = [int(a_id) for a_id in self.load_article_ids_from_journal(j_nlmid)]
            logging.debug('Paperboy.update_article_ids: downloaded={}'.format(
                article_ids))
            old_article_ids = self.cfg.journals_last_article_id_lists.get(
                j_nlmid, [])
            logging.debug('Paperboy.update_article_ids: old={}'.format(
                old_article_ids))
            new_article_ids = [int(a_id) for a_id in article_ids if a_id not in old_article_ids]
            logging.debug('Paperboy.update_article_ids: new={}'.format(
                new_article_ids))
            logging.info('Found {} new articles from {}.'.format(
                len(new_article_ids), j_abbr))

            # update config, storing all new articles since the last update
            self.cfg.journals_last_article_id_lists[j_nlmid] = new_article_ids

        # use today since only the date is stored, time not included, rely on
        # saved IDs to not show articles twice
        self.cfg.last_check_date = datetime.date.today()
        self.cfg.save_to_file()

    # load all articles
    def load_articles(self, article_ids):
        if len(article_ids) == 0:
            return

        parser = lxml.etree.XMLParser(remove_blank_text=True)

        id_string = ",".join(str(a_id) for a_id in article_ids)
        record = self.entrez_efetch(db='pubmed', retmode='xml', id=id_string)
        xml_tree = lxml.etree.fromstring(record, parser)

        for child in xml_tree:
            logging.debug('Paperboy.load_articles: {}'.format(
                lxml.etree.tostring(child)))
            if child.tag == 'PubmedArticle':
                self.articles.append(Article(child))
            else:
                logging.debug('Ignoring child node that is no PubmedArticle')

    # load the article data for all journals
    def load_all_articles(self):
        self.articles = []
        for j in self.cfg.journals:
            self.load_articles(self.cfg.journals_last_article_id_lists[j])

    # print the article data
    def show_articles(self, show_all):
        if show_all is True:
            articles_to_show = self.articles
            letters_and_editorials = 0
        else:
            pub_types_to_remove = ['Letter', 'Editorial', 'News', 'Interview']
            articles_to_show = []
            for a in self.articles:
                if len(set(pub_types_to_remove) & set(a.pub_type_list)) == 0:
                    articles_to_show.append(a)

            letters_and_editorials = len(self.articles) - len(articles_to_show)

        if len(articles_to_show) == 0:
            logging.info('No articles to show.')
        else:
            logging.info('Displaying {} articles{}.'.format(
                len(articles_to_show),
                '' if letters_and_editorials == 0 else
                ', excluding {} {}'.format(
                    letters_and_editorials, ', '.join(pub_types_to_remove)
                )
            ))

        for a in articles_to_show:
            print(a)

    def update_journal_list(self):
        if len(self.cfg.pubmed_journals) > 0 and self.cfg.pubmed_journals_last_update >= datetime.date.today() - datetime.timedelta(days=self.cfg.pubmed_journal_update_interval_days):
            return

        logging.info('Updating journal information.')
        logging.debug('Paperboy.update_journal_list: {}'.format(
            self.cfg.pubmed_journals_url))
        journals = []
        journal_divider_string = '--------------------------------------------------------'
        with urllib.request.urlopen(self.cfg.pubmed_journals_url) as response:
            content = response.read().decode('utf-8').splitlines()

            journal_content = []
            for line in content:
                if line.startswith(journal_divider_string[0:5]):
                    if len(journal_content) > 0:
                        journals.append(Journal(journal_content))
                        journal_content = []
                else:
                    journal_content.append(line)

        # save
        if len(journals) > 0:
            self.cfg.pubmed_journals = journals
            self.cfg.pubmed_journals_last_update = datetime.date.today()
            self.cfg.save_to_file()

    def get_journal(self, journal_nlmid):
        for pubmed_j in self.cfg.pubmed_journals:
            if journal_nlmid == pubmed_j.journal_data_dict[Journal.j_nlmid_key]:
                return pubmed_j
        return None

    def add_journal(self, journal_nlmid):
        self.update_journal_list()

        # get journal abbreviation for journal_nlmid
        j = self.get_journal(journal_nlmid)
        if j is None:
            logging.info('No journal with NlmId {} found.'.format(
                journal_nlmid))
            return

        j_abbr = j.journal_data_dict[Journal.j_medabbr_key]

        # check if journal_nlmid is already in journal list
        if journal_nlmid in self.cfg.journals:
            logging.info('{} (NlmId: {}) is already active.'.format(
                j_abbr, journal_nlmid))
            return

        # add journal_nlmid to journal list
        self.cfg.journals.append(journal_nlmid)
        self.cfg.save_to_file()
        logging.info('{} (NlmId: {}) was marked active.'.format(
            j_abbr, journal_nlmid))

    def remove_journal(self, journal_nlmid):
        self.update_journal_list()

        # get journal data for journal_nlmid
        j = self.get_journal(journal_nlmid)
        if j is None:
            logging.info('No journal with NlmId {} found.'.format(
                journal_nlmid))
            return

        j_abbr = j.journal_data_dict[Journal.j_medabbr_key]

        # check if journal_nlmid is in journal list
        if journal_nlmid in self.cfg.journals:
            logging.info('{} (NlmId: {}) is marked inactive.'.format(
                j_abbr, journal_nlmid))
            self.cfg.journals.remove(journal_nlmid)
            self.cfg.journals_last_article_id_lists.pop(journal_nlmid, None)
            self.cfg.save_to_file()
            return

        # already not on journal list, just note to user
        logging.info('{} (NlmId: {}) was already inactive.'.format(
            j_abbr, journal_nlmid))

    def list_active_journals(self):
        self.update_journal_list()

        # iterate only once over all pubmed journals
        active_journals = []
        for pubmed_j in self.cfg.pubmed_journals:
            for active_j in self.cfg.journals:
                if active_j == pubmed_j.journal_data_dict[Journal.j_nlmid_key]:
                    active_journals.append(pubmed_j)

        logging.info('Found {} active journals.'.format(len(active_journals)))
        for j in active_journals:
            print(j)

    def list_all_journals(self):
        self.update_journal_list()

        logging.info('Found {} journals.'.format(
            len(self.cfg.pubmed_journals)))
        for j in self.cfg.pubmed_journals:
            print('{}'.format(j))
