# /////////////////////////////////////////////////////////////// #
# !python3.6
# -*- coding: utf-8 -*-
# Python Script initially created on 2019-04-17 
# Compiled by Aly @ Grasselli Geomechanics Group, UofT, 2019 
# Created using PyCharm 
# Current Version - Dated Apr 23, 2018
# /////////////////////////////////////////////////////////////// #


# /// TIMER FUNCTION /// #

'''
Function to calculate the time
    Inputs:
        - Difference in time in seconds
    Returns:
        - Time in minutes and seconds
'''

def calc_timer_values(end_time):
    minutes, sec = divmod(end_time, 60)
    if end_time < 60:
        return ("\033[1m%.2f seconds\033[0m" % end_time)
    else:
        return ("\033[1m%d minutes and %d seconds\033[0m." % (minutes, sec))


# /// ADMINISTRATIVE AND SORTING OF FILES IN FOLDER /// #

''' 
Function to make a string bold for printing (on terminal)
    Inputs:
        - A string
    Returns:
        - Red / Green /Bold string
'''

def red_text(val):  # RED Bold text
    tex = "\033[1;31m%s\033[0m" % val
    return tex


def green_text(val):  # GREEN Bold text
    tex = "\033[1;92m%s\033[0m" % val
    return tex


def bold_text(val):  # Bold text
    tex = "\033[1m%s\033[0m" % val
    return tex


if __name__ == "__main__":
    try:
        red_text(text)
        green_text(text)
        bold_text(text)
    except KeyboardInterrupt:
        # print("\n\033[1;31;0mTERMINATED BY USER\n")
        exit("TERMINATED BY USER")

