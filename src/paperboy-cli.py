#!/usr/bin/env python3

import paperboy as pb
import sys
import logging
import argparse
import datetime


def func_up(args):
    p = pb.Paperboy()
    logging.info('Loading new articles since {}.'.format(
        p.cfg.last_check_date.strftime(pb.DATEFORMAT_USER)))
    p.update_article_ids()


def func_up_show(args):
    p = pb.Paperboy()
    logging.info('Loading new articles since {}.'.format(
        p.cfg.last_check_date.strftime(pb.DATEFORMAT_USER)))
    p.update_article_ids()
    p.load_all_articles()
    p.show_articles(args.show_all)


def func_show(args):
    p = pb.Paperboy()
    p.load_all_articles()
    p.show_articles(args.show_all)


def func_show_last_update(args):
    logging.info('Last update was on {}.'.format(
        pb.Paperboy().cfg.last_check_date.strftime(pb.DATEFORMAT_USER)))


def func_set_last_update(args):
    p = pb.Paperboy()
    if args.date[0].date() == p.cfg.last_check_date: 
        logging.info('New date {} is same as last update date; no changes needed.'.format(args.date[0].strftime(pb.DATEFORMAT_USER)))
        return

    logging.info('Setting last update date to {}.'.format(args.date[0].strftime(pb.DATEFORMAT_USER)))
    p.cfg.last_check_date = args.date[0].date()
    p.cfg.old_article_ids = []
    p.cfg.save_to_file();


def func_j_list_all(args):
    pb.Paperboy().list_all_journals()


def func_j_list_active(args):
    pb.Paperboy().list_active_journals()


def func_j_add(args):
    pb.Paperboy().add_journal(args.nlmid[0])


def func_j_remove(args):
    pb.Paperboy().remove_journal(args.nlmid[0])


def valid_date(date_str):
    try:
        return datetime.datetime.strptime(date_str, pb.DATEFORMAT_USER)
    except ValueError:
        msg = 'Invalid date: {}'.format(date_str)
        raise argparse.ArgumentTypeError(msg)


def main():
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
