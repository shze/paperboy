#!/usr/bin/env python3

from PyQt5.QtWidgets import QApplication, QSystemTrayIcon, QMenu, QAction
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import QTimer, QTime, QThread
from notifypy import Notify
import paperboy as pb
import logging
import time
import os


MENU_SHOW_TEXT_NONE = "No new articles"
MENU_SHOW_TEXT_NEW = "Show {} new articles"
MENU_SHOW_TEXT_NEW_RECENT = "Show {} new and {} recent articles"
MENU_SHOW_TEXT_RECENT = "Show {} recent articles"

MENU_UPDATE_TEXT_IDLE = "Check for new articles"
MENU_UPDATE_TEXT_UPDATING = "Checking for new articles.."


# pyxdg does not have the xdg_state_home (yet); use this for now
def get_xdg_state_home():
    env_xdg_state_home = os.environ.get('XDG_STATE_HOME')
    env_home = os.path.expanduser("~")
    if env_xdg_state_home is None or not env_xdg_state_home.strip():
        return os.path.join(env_home, ".local", "state")
    else:
        return env_xdg_state_home.strip()


class UpdateWorker(QThread):
    def run(self): 
        p = pb.Paperboy()
        p.update_article_ids()

        notification = Notify(
            default_notification_title="Function Message",
            #default_application_name="Great Application",
            default_notification_application_name="not a python app",
            default_notification_icon="/usr/share/icons/hicolor/32x32/status/laptoptrusted.svg",
            #default_notification_audio="path/to/sound.wav"
        )
        notification.message = "Even cooler message."
        notification.send()

        time.sleep(10)

        self.finished.emit()


class App:
    def __init__(self):
        self.init_logging()

        self.app = QApplication([])
        self.app.setQuitOnLastWindowClosed(False)
        self.init_tray_icon()

        self.timer = QTimer()
        self.timer.timeout.connect(self.f_update_articles)
        self.timer.start(1000)

    def init_logging(self):
        # set up logging
        log_path = get_xdg_state_home()
        logfile = os.path.join(log_path, pb.PROJECTNAME + ".log")
        # make sure dirs exist
        os.makedirs(log_path, exist_ok=True)

        # always log to terminal: only INFO level and short format
        h_term = logging.StreamHandler()
        h_term.setLevel(logging.INFO)
        h_term.setFormatter(logging.Formatter(fmt='%(levelname)s: %(message)s'))

        debug = False # TODO needs to be in args
        if debug:
            # log to file: DEBUG level and long format
            h_file = logging.FileHandler(logfile)
            h_file.setLevel(logging.DEBUG)
            h_file.setFormatter(logging.Formatter(fmt='%(asctime)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))

            logging.basicConfig(level=logging.DEBUG, handlers=[h_term, h_file])
            logging.info('Logging debug info to {}.'.format(logfile))
        else:
            logging.basicConfig(level=logging.INFO, handlers=[h_term])

    def init_tray_icon(self):
        icon = QIcon.fromTheme("applications-science")
        self.tray = QSystemTrayIcon()
        self.tray.setIcon(icon)
        self.tray.setVisible(True)

        self.menu = QMenu()

        self.m_show_articles = QAction(MENU_SHOW_TEXT_NONE)
        # 3 states: a) disabled, no articles. b) enabled, old articles. c) enabled, new articles.
        self.m_show_articles.setEnabled(False)
        self.m_show_articles.setIcon(QIcon.fromTheme('text-plain'))
        self.menu.addAction(self.m_show_articles)

        self.m_update_articles = QAction(MENU_UPDATE_TEXT_IDLE)
        self.m_update_articles.triggered.connect(self.f_update_articles)
        self.m_update_articles.setIcon(QIcon.fromTheme('view-refresh'))
        self.menu.addAction(self.m_update_articles)

        self.m_quit = QAction("Quit")
        self.m_quit.triggered.connect(self.app.quit)
        self.m_quit.setIcon(QIcon.fromTheme('application-exit'))
        self.menu.addAction(self.m_quit)

        self.tray.setContextMenu(self.menu)

    def run(self):
        self.app.exec_()

    def f_update_running(self, is_running):
        if is_running:
            self.timer.stop()
            self.m_update_articles.setText(MENU_UPDATE_TEXT_UPDATING)
            self.m_update_articles.setEnabled(False)
        else:
            self.timer.start(600000) # 10min
            self.m_update_articles.setText(MENU_UPDATE_TEXT_IDLE)
            self.m_update_articles.setEnabled(True)

    def f_update_articles(self):
        self.f_update_running(True)

        self.w = UpdateWorker()
        self.w.finished.connect(lambda: self.f_update_running(False))
        self.w.start()


if __name__ == '__main__':
    App().run()
