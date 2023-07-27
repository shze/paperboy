#!/usr/bin/env python3

from PyQt5.QtWidgets import QApplication, QSystemTrayIcon, QMenu, QAction
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import QTimer, QTime
from notifypy import Notify
import paperboy as pb
import logging

def f_update_articles():
    global timer
    timer.stop()
    # TODO disable check for new articles menu
    
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

    # TODO enable menu again
    timer.start(600000) # 10min

# set up logging
debug = False # TODO needs to be in args
logging.basicConfig(format='%(levelname)s: %(message)s',
                    level=logging.DEBUG if debug is True else logging.INFO)

app = QApplication([])
app.setQuitOnLastWindowClosed(False)

# TODO better icon, depending on DE and state
icon = QIcon.fromTheme("applications-science")
tray = QSystemTrayIcon()
tray.setIcon(icon)
tray.setVisible(True)

menu = QMenu()

show_articles = QAction("No new articles")
# 4 states: a) disabled, no articles. b) enabled, old articles. c) enabled, new articles. d) checking.
show_articles.setEnabled(False)
show_articles.setIcon(QIcon.fromTheme('system-software-install'))
menu.addAction(show_articles)

update_articles = QAction("Check for new articles")
update_articles.triggered.connect(f_update_articles)
update_articles.setIcon(QIcon.fromTheme('view-refresh'))
menu.addAction(update_articles)

quit = QAction("Quit")
quit.triggered.connect(app.quit)
quit.setIcon(QIcon.fromTheme('application-exit'))
menu.addAction(quit)

tray.setContextMenu(menu)

timer = QTimer()
timer.timeout.connect(f_update_articles)
timer.start(1000)

app.exec_()
