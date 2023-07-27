#!/usr/bin/env python3

from PyQt5.QtWidgets import QApplication, QSystemTrayIcon, QMenu, QAction
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import QTimer, QTime, QThread
from notifypy import Notify
import paperboy as pb
import logging
import time


MENU_UPDATE_TEXT_IDLE = "Check for new articles"
MENU_UPDATE_TEXT_UPDATING = "Checking for new articles.."


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
        # set up logging
        debug = False # TODO needs to be in args
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.DEBUG if debug is True else logging.INFO)

        self.app = QApplication([])
        self.app.setQuitOnLastWindowClosed(False)

        icon = QIcon.fromTheme("applications-science")
        self.tray = QSystemTrayIcon()
        self.tray.setIcon(icon)
        self.tray.setVisible(True)

        self.menu = QMenu()

        self.m_show_articles = QAction("No new articles")
        # 4 states: a) disabled, no articles. b) enabled, old articles. c) enabled, new articles. d) checking.
        self.m_show_articles.setEnabled(False)
        self.m_show_articles.setIcon(QIcon.fromTheme('system-software-install'))
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

        self.timer = QTimer()
        self.timer.timeout.connect(self.f_update_articles)
        self.timer.start(1000)

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
