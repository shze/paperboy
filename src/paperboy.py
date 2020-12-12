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

APPNAME = 'paperboy'
CONFFILE = '{}.conf'.format(APPNAME)
APPAUTHOR = 'heinze'
EMAIL = 'sten.heinze@gmail.de'

# https://stackoverflow.com/questions/11875770/how-to-overcome-datetime-datetime-not-json-serializable
def datetime_serializer(obj):
    if isinstance(obj, (datetime, date)):
        return obj.strftime('%Y-%m-%d')
    raise TypeError('Type {} not JSON serializable'.format(type(obj)))

# https://stackoverflow.com/questions/8793448/how-to-convert-to-a-python-datetime-object-with-json-loads
def datetime_parser(json_dict):
    for (key, value) in json_dict.items():
        try:
            json_dict[key] = datetime.strptime(value, '%Y-%m-%d')
        except (ValueError, TypeError):
            pass
    return json_dict

class Config:
    def __init__(self):
        self.last_check_date = date.today() - timedelta(days = 4)
        self.email = EMAIL
        #self.pubmed_journals_url = "https://ftp.ncbi.nih.gov/pubmed/J_Medline.txt"
        self.journals = ['science']
        self.journals_last_article_id_lists = dict() # dict {journal_string : [article_id, ..], ..}
        
        self.load_from_file()
        
    def get_member_dict(self):
        member_dict = {attr : getattr(self, attr) for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__")}
        logging.debug('Config.get_member_dict: {}'.format(member_dict))
        return member_dict
    
    def set_members(self, member_dict):
        members = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__")]
        for m in members:
            if m in member_dict:
                setattr(self, m, member_dict[m])
        logging.debug('Config.set_members: {}'.format(member_dict))
        
    def load_from_file(self):
        cfg_file = os.path.join(appdirs.user_config_dir(APPNAME), CONFFILE)
        
        try:
            with open(cfg_file, 'r') as f:
                self.set_members(json.load(f, object_hook = datetime_parser))
                logging.debug('Config.load_from_file: loading from {}'.format(cfg_file))
        except:
            logging.debug('Config.load_from_file: could not read {}, using defaults'.format(cfg_file))
            raise
        
    def save_to_file(self):
        folder = appdirs.user_config_dir(APPNAME)
        os.makedirs(folder, exist_ok = True)
        
        cfg_file = os.path.join(folder, CONFFILE)
        logging.debug('Config.save_to_file: saving to {}'.format(cfg_file))
        
        with open(cfg_file, 'w') as f:
            json.dump(self.get_member_dict(), f, default = datetime_serializer, indent=4, sort_keys=True)
  
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
        #xpath_author_firstname = './Author/ForeName'
        #xpath_author_middle = './Author/Initials'
        xpath_author_collective_name = './Author/CollectiveName'
        
        self.missing_data = []
        self.doi = 'https://doi.org/{}'.format(self.extract_first_data(article_xml_tree, 'DOI', xpath_doi))
        self.pmid = self.extract_first_data(article_xml_tree, 'PMID', xpath_pmid)
        self.title = self.extract_first_data(article_xml_tree, 'Title', xpath_title)
        self.pub_type_list = self.extract_all_data_list(article_xml_tree, 'PubType', xpath_pub_type)
        self.journal_abbrev = self.extract_first_data(article_xml_tree, 'JournalAbbrev', xpath_journal_abbrev)
        self.journal_vol = self.extract_first_data(article_xml_tree, 'JournalVol', xpath_journal_vol)
        self.journal_issue = self.extract_first_data(article_xml_tree, 'JournalIssue', xpath_journal_issue)
        
        y = self.extract_first_data(article_xml_tree, 'JournalYear', xpath_journal_y)
        m = self.extract_first_data(article_xml_tree, 'JournalMonth', xpath_journal_m)
        d = self.extract_first_data(article_xml_tree, 'JournalDay', xpath_journal_d)
        logging.debug('Article.__init__: y={} m={} d={}'.format(y, m, d))
        # some journals do not use day or month field, set it 1 to avoid breaking the date parsing
        if m == None:
            m = 1
        if d == None:
            d = 1
        try:
            self.journal_date = datetime.strptime('{}-{}-{}'.format(y, m, d), '%Y-%b-%d')
        except:
            try:
                self.journal_date = datetime.strptime('{}-{}-{}'.format(y, m, d), '%Y-%m-%d')
            except:
                self.missing_data.append('JournalDate')
        
        authorlist_xml_tree_list = article_xml_tree.xpath(xpath_authorlist)
        if not len(authorlist_xml_tree_list) == 1:
            self.missing_data.append('AuthorList')
            self.author_lastname_list = []
        else:
            authorlist_xml_tree = authorlist_xml_tree_list[0]
            self.author_lastname_list = self.extract_all_data_list(authorlist_xml_tree, 'Lastname', xpath_author_lastname)
            if self.author_lastname_list == None:
                self.author_lastname_list = self.extract_all_data_list(authorlist_xml_tree, 'CollectiveName', xpath_author_collective_name)
        logging.debug('Article.__init__: authors={}'.format(self.author_lastname_list))
        
        logging.debug('Article.__init__: {}'.format(self.__repr__()))
        if len(self.missing_data) > 0:
            logging.debug('Missing data for PMID {}: {}'.format(self.pmid, ', '.join(self.missing_data)))

    # object representation
    def __repr__(self):
        return 'Article(DOI={}, PMID={}, \'{}\' {}. {}, {}, {}, {}, {})' \
            .format(self.doi, self.pmid, self.title, ', '.join(self.author_lastname_list), \
                    self.journal_abbrev, self.journal_vol, self.journal_issue, \
                    self.journal_date.strftime('%Y-%m-%d'), ', '.join(self.pub_type_list))
    
    # user readable object information that fits on one terminal line, min 110 chars
    def __str__(self):
        # use 110 as minimum width, b/c the title would be shortened to zero chars if we're trying to fit everything in 80 chars
        max_len = max(110, os.get_terminal_size()[0])
        
        # some article have no authors, e.g. news, but are still categorized as 'Journal Article'
        author_lastnames = ''
        if len(self.author_lastname_list) > 0:
            author_lastnames = '{}. '.format(', '.join(self.author_lastname_list))
            
        # some article have many publication types, remove less informative ones
        pub_types_to_remove = ['Research Support, Non-U.S. Gov\'t', 
                               'Research Support, U.S. Gov\'t, Non-P.H.S.', 
                               'Research Support, U.S. Gov\'t, P.H.S.',
                               'Research Support, N.I.H., Extramural',
                               'Research Support, N.I.H., Intramural']
        pub_type_list_to_show = [p for p in self.pub_type_list if not p in pub_types_to_remove]
            
        journal_doi_pmid = '{} ('.format(self.journal_abbrev) \
            + Fore.YELLOW + '{}'.format(self.doi) \
            + Style.RESET_ALL + ') (PMID ' \
            + Fore.YELLOW + '{}'.format(self.pmid) \
            + Style.RESET_ALL + ') ({})'.format(', '.join(pub_type_list_to_show))
        
        title = Style.BRIGHT + '{}'.format(self.title) + Style.RESET_ALL + ' '
        
        result = title + author_lastnames + journal_doi_pmid
        
        # https://stackoverflow.com/questions/2186919/getting-correct-string-length-in-python-for-strings-with-ansi-color-codes
        ESC = Literal('\x1b')
        integer = Word(nums)
        escapeSeq = Combine(ESC + '[' + Optional(delimitedList(integer, ';')) + oneOf(list(alphas)))
        nonAnsiString = lambda s: Suppress(escapeSeq).transformString(s)

        if len(nonAnsiString(result)) > max_len and len(self.author_lastname_list) > 2:
            logging.debug('Article.__str__: len_str={} > max_len={}, len_auth={} > 2'.format(len(nonAnsiString(result)), max_len, len(self.author_lastname_list)))
            author_lastnames = ', '.join(self.author_lastname_list[0:1]) + ', ({} more). '.format(len(self.author_lastname_list) - 1)
            result = title + author_lastnames + journal_doi_pmid
                
        if len(nonAnsiString(result)) > max_len:
            logging.debug('Article.__str__: len_str={} > max_len={}'.format(len(nonAnsiString(result)), max_len))
            title_len = len(self.title) + max_len - len(nonAnsiString(result)) - 2
            title = Style.BRIGHT + '{}..'.format(self.title[0:title_len].strip()) + Style.RESET_ALL + ' '
            result = title + author_lastnames + journal_doi_pmid
        
        return result
        
    def extract_first_data(self, xml_tree, description, xpath_spec):
        data_list = xml_tree.xpath(xpath_spec)
        if len(data_list) == 1:
            logging.debug('Article.extract_first_data: {}'.format(etree.tostring(data_list[0])))
            
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
            logging.debug('Article.extract_all_data_list: {}'.format([etree.tostring(d) for d in data_list]))
            
            data_string_list = [''.join(d.itertext()) for d in data_list]
            logging.debug('Article.extract_all_data_list: {}'.format(data_string_list))
            return [d.replace('\r', '').replace('\n', '') for d in data_string_list]
        else:
            self.missing_data.append(description)
            return None

class Paperboy:
    # constructor
    def __init__(self):
        self.cfg = Config()
        self.articles = []
        Entrez.email = self.cfg.email
        
    def entrez_esearch(self, *args, **kwargs):
        try:
            handle = Entrez.esearch(*args, **kwargs)
            record = Entrez.read(handle)
            handle.close()
            return record 
        except (IOError, OSError):
            logging.error('Paperboy.entrez_esearch: I/O Error searching data from NCBI Entrez.')
            sys.exit()
        
    def entrez_efetch(self, *args, **kwargs):
        try:
            handle = Entrez.efetch(*args, **kwargs)
            record = handle.read()
            handle.close()
            return record 
        except (IOError, OSError):
            logging.error('Paperboy.entrez_efetch: I/O Error fetching data from NCBI Entrez.')
            sys.exit()
        
    # load and return all (up to 100000) article IDs for a given journal
    def load_article_ids_from_journal(self, journal):
        term = journal + '[ta] AND ' \
                + self.cfg.last_check_date.strftime('%Y/%m/%d') + ":" \
                + date.today().strftime('%Y/%m/%d') + '[dp]'
        record = self.entrez_esearch(db = 'pubmed', term = term)
        logging.debug('Paperboy.load_article_ids_from_journal: {}'.format(record))
            
        article_ids = record['IdList']
            
        # limit to 100000, otherwise needs retstart loop
        # https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
        count = min(int(record['Count']), 100000) 
        retsize = len(article_ids)
            
        if count > retsize:
            logging.debug('Paperboy.load_article_ids_from_journal: count={} > retsize={}'.format(count, retsize))
            record = self.entrez_esearch(db = 'pubmed', term = term, retmax = count, retstart = retsize)
            logging.debug('Paperboy.load_article_ids_from_journal: {}'.format(record))
                
            article_ids.extend(record['IdList'])
            
        return article_ids
    
    # load all articles
    def load_articles(self, article_ids):
        if len(article_ids) == 0:
            return
        
        parser = etree.XMLParser(remove_blank_text = True)
        
        id_string = ",".join(str(a_id) for a_id in article_ids)
        record = self.entrez_efetch(db = 'pubmed', retmode = 'xml', id = id_string)
        xml_tree = etree.fromstring(record, parser)
        
        for child in xml_tree:
            logging.debug('Paperboy.load_articles: {}'.format(etree.tostring(child)))
            if child.tag == 'PubmedArticle':
                self.articles.append(Article(child))
            else:
                logging.debug('Ignoring child node that is not a PubmedArticle')

    # load article IDs for all journals
    def update_article_ids(self):
        logging.info('Loading new articles since {}.'.format(self.cfg.last_check_date.strftime('%Y/%m/%d')))
        
        for j in self.cfg.journals:
            # load article IDs
            article_id_strings = self.load_article_ids_from_journal(j)
            logging.debug('Paperboy.update_article_ids: {}'.format(article_id_strings))
            new_article_ids = [int(a_id) for a_id in article_id_strings if a_id not in self.cfg.journals_last_article_id_lists.get(j, [])]
            logging.debug('Paperboy.update_article_ids: {}'.format(new_article_ids))
            logging.info('Found {} new articles from {}.'.format(len(new_article_ids), j))
            
            # update config
            self.cfg.journals_last_article_id_lists[j] = new_article_ids
            
        self.cfg.last_check_date = date.today()
        self.cfg.save_to_file()
        
    # load the article data for all journals
    def load_all_articles(self):
        self.articles = []
        for j in self.cfg.journals:
            self.load_articles(self.cfg.journals_last_article_id_lists[j])
            
    ## print the article data 
    def show_articles(self, show_all):
        if show_all == True:
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
            logging.info('Displaying {} articles{}.'.format( \
                len(articles_to_show), \
                '' if letters_and_editorials == 0 else ', excluding {} {}'.format(letters_and_editorials, ', '.join(pub_types_to_remove)) \
            ))
            
        for a in articles_to_show:
            print(a)
        
    def add_journal(self):
        pass
    
    def remove_journal(self):
        pass
    
    def list_journals(self):
        pass

def main():
    init() # for colorama
    
    # cli
    cmd_up = 'update'
    cmd_up_show = 'update-show'
    cmd_show = 'show'
    parser = argparse.ArgumentParser()
    parser.add_argument('command',
        help = 'command to be executed (default: %(default)s)',
        choices = [cmd_up, cmd_up_show, cmd_show], 
        nargs = '?',
        default = cmd_up_show
    )
    parser.add_argument('-s', '--show-all', action = 'store_true', help = 'show also letters and editorials')
    parser.add_argument('-d', '--debug', action = 'store_true', help = 'show debug output')
    args = parser.parse_args()
    
    # set up logging
    logging.basicConfig(format = '%(levelname)s: %(message)s', level = logging.DEBUG if args.debug == True else logging.INFO) 
    
    p = Paperboy()
    
    if args.command == cmd_up:
        p.update_article_ids()
    elif args.command == cmd_up_show:
        p.update_article_ids()
        p.load_all_articles()
        p.show_articles(args.show_all)
    elif args.command == cmd_show:
        p.load_all_articles()
        p.show_articles(args.show_all)

if __name__ == '__main__':
    main()

# TODO EN locale for parsing pubmed date https://stackoverflow.com/questions/985505/locale-date-formatting-in-python
# TODO cli: show pmid, add, remove, list
# TODO rate limit
# TODO format+document source according to py rules
