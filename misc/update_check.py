# -*- encoding: utf-8 -*-
# möp
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
from __future__ import print_function

urlprefix = "https://www.xs3-data.uni-freiburg.de/structurefinder"

from urllib import urlopen


def get_current_strf_version(VERSION=0, silent=True):
    """
    determines the current version of DSR on the web server

    >>> get_current_strf_version()
    '202'

    Returns
    -------
    version number
    :type: str
    """
    import socket
    socket.setdefaulttimeout(3)
    try:
        response = urlopen('{}/version.txt'.format(urlprefix))
    except IOError:
        if not silent:
            print("*** Unable to connect to update server. No Update possible. ***")
        return 0
    try:
        version = response.readline().decode('ascii').strip()
    except ValueError:
        return 0
    return version


def is_update_needed(VERSION=0, silent=True):
    """
    Decides if an update of DSR is needed.
    :return: bool
    >>> is_update_needed()
    False
    """
    version = get_current_strf_version(VERSION=VERSION, silent=True)
    if int(VERSION) < int(version):
        return True
    else:
        return False