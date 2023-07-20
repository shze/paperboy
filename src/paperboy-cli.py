#!/usr/bin/env python

from datetime import date, datetime, timedelta
from Bio import Entrez
from lxml import etree
from colorama import init, Fore, Back, Style
from pyparsing import *
import os
import sys
import logging
import argparse
import appdirs
import json
import urllib.request
import re

APPNAME = 'paperboy'
CONFFILE = '{}.conf'.format(APPNAME)
APPAUTHOR = 'heinze'
EMAIL = 'sten.heinze@gmail.de'
DATEFORMAT_INTERNAL = '%Y-%m-%d'
DATEFORMAT_ENTREZ = '%Y/%m/%d'
DATEFORMAT_USER = '%Y-%m-%d'
# use min width so that output is not shortened so much to become unreadable
MIN_TERMINAL_WIDTH = 110
try:
    MAX_LEN = max(MIN_TERMINAL_WIDTH, os.get_terminal_size()[0])
except OSError:
    MAX_LEN = sys.maxsize

# https://stackoverflow.com/questions/2186919/getting-correct-string-length-in-python-for-strings-with-ansi-color-codes
ESC = Literal('\x1b')
integer = Word(nums)
escapeSeq = Combine(
    ESC + '[' + Optional(delimitedList(integer, ';')) + oneOf(list(alphas)))


def nonAnsiString(s): return Suppress(escapeSeq).transformString(s)


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
            return (date.fromisoformat(s), end)
        else:
            return (s, end)

    def object_hook_func(self, obj):
        if '__type__' in obj and obj['__type__'] == 'Journal':
            return Journal(obj['__data__'])
        return obj


def custom_encoder(obj):
    if isinstance(obj, date):
        return obj.isoformat()
    if isinstance(obj, Journal):
        return {'__type__': 'Journal', '__data__': obj.medline_format()}
    raise TypeError('Type {} not JSON serializable'.format(type(obj)))


class Config:
    def __init__(self):
        self.last_check_date = date.today() - timedelta(days=4)
        self.email = EMAIL
        self.pubmed_journals_url = "https://ftp.ncbi.nih.gov/pubmed/J_Medline.txt"
        self.pubmed_journals = []
        self.pubmed_journals_last_update = date.min
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
        cfg_file = os.path.join(appdirs.user_config_dir(APPNAME), CONFFILE)

        try:
            with open(cfg_file, 'r') as f:
                self.set_members(json.load(f, cls=JSONDecoder2))
                logging.debug(
                    'Config.load_from_file: loading from {}'.format(cfg_file))
        except OSError:
            logging.debug(
                'Config.load_from_file: could not read {}, using defaults'.format(cfg_file))

    def save_to_file(self):
        folder = appdirs.user_config_dir(APPNAME)
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
            self.journal_date = datetime.strptime(
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
            + Fore.YELLOW + '{}'.format(self.doi) + Style.RESET_ALL \
            + ') (PMID ' \
            + Fore.YELLOW + '{}'.format(self.pmid) + Style.RESET_ALL \
            + ') ({})'.format(', '.join(pub_type_list_to_show))

        title = Style.BRIGHT + '{}'.format(self.title) + Style.RESET_ALL + ' '

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
            title = Style.BRIGHT \
                + '{}..'.format(self.title[0:title_len].strip()) \
                + Style.RESET_ALL + ' '
            result = title + author_lastnames + journal_doi_pmid

        return result

    def extract_first_data(self, xml_tree, description, xpath_spec):
        data_list = xml_tree.xpath(xpath_spec)
        if len(data_list) == 1:
            logging.debug('Article.extract_first_data: {}'.format(
                etree.tostring(data_list[0])))

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
                [etree.tostring(d) for d in data_list]))
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
            + Fore.YELLOW \
            + '{}'.format(self.journal_data_dict[Journal.j_nlmid_key]) \
            + Style.RESET_ALL + ')' \
            if Journal.j_nlmid_key in self.journal_data_dict else ""

        issn_list = []
        if Journal.j_issn_print_key in self.journal_data_dict:
            issn_list.append(self.journal_data_dict[Journal.j_issn_print_key])
        if Journal.j_issn_online_key in self.journal_data_dict:
            issn_list.append(self.journal_data_dict[Journal.j_issn_online_key])
        issn_part = ""
        if len(issn_list) > 0:
            issn_part = ' (ISSN: ' + ', '.join(issn_list) + ')'

        result = 'Journal: ' + Style.BRIGHT + \
            '{}'.format(title) + Style.RESET_ALL + \
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
        Entrez.email = self.cfg.email

    def entrez_esearch(self, *args, **kwargs):
        logging.debug('Paperboy.entrez_esearch: {} {}'.format(args, kwargs))

        try:
            handle = Entrez.esearch(*args, **kwargs)
            record = Entrez.read(handle)
            handle.close()
            return record
        except (IOError, OSError):
            logging.error(
                'Paperboy.entrez_esearch: I/O Error searching data from NCBI Entrez.')
            sys.exit()

    def entrez_efetch(self, *args, **kwargs):
        logging.debug('Paperboy.entrez_efetch: {} {}'.format(args, kwargs))

        try:
            handle = Entrez.efetch(*args, **kwargs)
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
            + date.today().strftime(DATEFORMAT_ENTREZ) + '[dp]'
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
        logging.info('Loading new articles since {}.'.format(
            self.cfg.last_check_date.strftime(DATEFORMAT_USER)))

        active_journals = []
        for pubmed_j in self.cfg.pubmed_journals:
            for active_j in self.cfg.journals:
                if active_j == pubmed_j.journal_data_dict[Journal.j_nlmid_key]:
                    active_journals.append(pubmed_j)

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
        self.cfg.last_check_date = date.today()
        self.cfg.save_to_file()

    # load all articles
    def load_articles(self, article_ids):
        if len(article_ids) == 0:
            return

        parser = etree.XMLParser(remove_blank_text=True)

        id_string = ",".join(str(a_id) for a_id in article_ids)
        record = self.entrez_efetch(db='pubmed', retmode='xml', id=id_string)
        xml_tree = etree.fromstring(record, parser)

        for child in xml_tree:
            logging.debug('Paperboy.load_articles: {}'.format(
                etree.tostring(child)))
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
        if len(self.cfg.pubmed_journals) > 0 and self.cfg.pubmed_journals_last_update >= date.today() - timedelta(days=self.cfg.pubmed_journal_update_interval_days):
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
            self.cfg.pubmed_journals_last_update = date.today()
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


def func_up(args):
    Paperboy().update_article_ids()


def func_up_show(args):
    p = Paperboy()
    p.update_article_ids()
    p.load_all_articles()
    p.show_articles(args.show_all)


def func_show(args):
    p = Paperboy()
    p.load_all_articles()
    p.show_articles(args.show_all)


def func_show_last_update(args):
    logging.info('Last update was on {}.'.format(
        Paperboy().cfg.last_check_date.strftime(DATEFORMAT_USER)))


def func_set_last_update(args):
    p = Paperboy()
    if args.date[0].date() == p.cfg.last_check_date: 
        logging.info('New date {} is same as last update date; no changes needed.'.format(args.date[0].strftime(DATEFORMAT_USER)))
        return

    logging.info('Setting last update date to {}.'.format(args.date[0].strftime(DATEFORMAT_USER)))
    p.cfg.last_check_date = args.date[0].date()
    p.cfg.old_article_ids = []
    p.cfg.save_to_file();


def func_j_list_all(args):
    Paperboy().list_all_journals()


def func_j_list_active(args):
    Paperboy().list_active_journals()


def func_j_add(args):
    Paperboy().add_journal(args.nlmid[0])


def func_j_remove(args):
    Paperboy().remove_journal(args.nlmid[0])


def valid_date(date_str):
    try:
        return datetime.strptime(date_str, DATEFORMAT_USER)
    except ValueError:
        msg = 'Invalid date: {}'.format(date_str)
        raise argparse.ArgumentTypeError(msg)


def main():
    init()  # for colorama

    # cli
    cmd_up_show = 'update-show'
    cmd_show = 'show'
    cmd_set_last_up = 'set-last-update'
    cmd_j_add = 'journal-add'
    cmd_j_remove = 'journal-remove'
    cmd_dict = {'update': func_up, 
                cmd_up_show: func_up_show,
                cmd_show: func_show,
                'show-last-update': func_show_last_update,
                cmd_set_last_up: func_set_last_update,
                'journal-list-all': func_j_list_all,
                'journal-list-active': func_j_list_active,
                cmd_j_add: func_j_add,
                cmd_j_remove: func_j_remove}

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', action='store_true',
                        help='show debug output')
    subparsers = parser.add_subparsers(
        help='command to be executed (default: {})'.format(cmd_up_show),
        dest='selected_cmd')
    cmd_parser = {}
    for cmd, func in cmd_dict.items():
        cmd_parser[cmd] = subparsers.add_parser(cmd)
        cmd_parser[cmd].set_defaults(func=func)

    cmd_parser[cmd_up_show].add_argument(
        '-s', '--show-all', action='store_true',
        help='show also letters and editorials')
    cmd_parser[cmd_show].add_argument(
        '-s', '--show-all', action='store_true',
        help='show also letters and editorials')
    cmd_parser[cmd_set_last_up].add_argument('date', nargs=1, type=valid_date)
    cmd_parser[cmd_j_add].add_argument('nlmid', nargs=1, type=str)
    cmd_parser[cmd_j_remove].add_argument('nlmid', nargs=1, type=str)

    # if no command was passed, append default
    if len(set(sys.argv) & set(cmd_dict.keys())) == 0:
        sys.argv.append(cmd_up_show)

    args = parser.parse_args()

    # set up logging
    logging.basicConfig(format='%(levelname)s: %(message)s',
                        level=logging.DEBUG if args.debug is True else logging.INFO)
    # call function for command
    args.func(args)


if __name__ == '__main__':
    main()

# TODO rate limit
# TODO document source according to py rules
