#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Function definitions common to all programs.


"""

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## libraries to use

import re
import os
import logging
import time

import getpass
import logging

from subprocess import Popen, PIPE
from email.mime.text import MIMEText



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function definitions


'''
creates a logging instance
'''
def get_logger(log_name, log_file, log_level = "INFO"):

    log = logging.getLogger(log_name)
    handler = logging.FileHandler(log_file)
    formatter = logging.Formatter('%(filename)-15s:%(process)d %(asctime)s %(levelname)s: %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(log_level)

    return log

'''
Checkpoint the status plus a timestamp
- appends the status

@param status_log: /path/to/status.log (or whatever you name it)
@param status: status to append to status.log

'''
def checkpoint_step(status_log, status):

    status_line = "%s,%s\n" % (status, time.strftime("%Y-%m-%d %H:%M:%S"))

    with open(status_log, "a") as myfile:
        myfile.write(status_line)

'''
returns the last step (status) from the pipeline
@param status_log: /path/to/status.log (or whatever you name it)
@param log: logger object

@return last status in the status log, "start" if nothing there
'''
def get_status(status_log, log = None):

    #status_log = "%s/%s" % (output_path, "test_status.log")

    status = "start"
    timestamp = str(time.strftime("%Y-%m-%d %H:%M:%S"))

    if os.path.isfile(status_log):

        fh = open(status_log, 'r')
        lines = fh.readlines()
        fh.close()

        for line in lines:
            line_list = line.split(",")
            assert(len(line_list) == 2)
            status = str(line_list[0]).strip()
            timestamp = str(line_list[1]).strip()

        if not status:
            status = "start"

        if log:
            log.info("Last checkpointed step: %s (%s)", status, timestamp)

    else:

        if log:
            log.info("Cannot find status.log (%s), assuming new run", status_log)

    status = status.strip().lower()

    return status



'''
run a command from python

@param cmd: command to run
@param live: False = run in dry mode (print command), True = run normally
@param log: logger object

@return std_out, std_err, exit_code
'''
def run_command(cmd, live = False, log = None):

    # our outputs
    std_out = None
    std_err = None
    exit_code = None

    if log:
        log.info("Running %s", cmd)

    if cmd:

        if not live:
            std_out = "Not live: cmd = '%s'" % (cmd)
            exit_code = 0
        else:

            p = Popen(cmd, stdout = PIPE, stderr = PIPE, shell = True)
            #p.wait()
            std_out, std_err = p.communicate()
            exit_code = p.returncode

    # use post_mortem_cmd(cmd, exit_code, std_out, std_err, log) to print this
    #if log:
        #log.info("Return vals from " + str(cmd) + " are: exitCode " + str(exit_code) + " stdOut " + str(std_out) + " stdErr " + str(std_err))


    return std_out, std_err, exit_code

'''
returns human readable file size
@param num = file size (e.g. 1000)

@return: readable float e.g. 1.5 KB
'''
def human_size(num):

    for x in ['bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'XB']:
        if (num < 1024.0):
            return "%3.1f %s" % (num, x)
        num /= 1024.0

    return "%3.1f %s" % (num, 'ZB')



'''
send out email
@param emailTo: email receipient (e.g. bryce@lbl.gov)
@param emailSubject: subject line for the email
@param emailBody: content of the email
@param emailFrom: optional email from
- Alex Boyd has a new version of this
'''
def send_email(emailTo, emailSubject, emailBody, emailFrom = 'rqc@jgi-psf.org', log = None):

    msg = ""
    errFlag = 0

    if not emailTo:
        msg = "- send_email: email_to parameter missing!"

    if not emailSubject:
        msg = "- send_email: email_subject parameter missing!"

    if not emailBody:
        msg = "- send_email: email_body parameter missing!"


    if errFlag == 0:
        msg = "- sending email to: %s" % (emailTo)

    if log:
        log.info(msg)
    else:
        print msg

    if errFlag == 1:
        return 0

    # assume html
    emailMsg = MIMEText(emailBody, "html") # vs "plain"
    emailMsg['Subject'] = emailSubject
    emailMsg['From'] = emailFrom
    emailMsg['To'] = emailTo

    p = Popen(["/usr/sbin/sendmail", "-t"], stdin = PIPE)
    p.communicate(emailMsg.as_string())

    return errFlag

'''
Write to rqc_file (e.g. rqc-files.tmp) the file_key and file_value

@param rqc_file_log: full path to file containing key=file
@param file_key: key for the entry
@param file_value: value for the entry

'''
def append_rqc_file(rqc_file_log, file_key, file_value, log=None):

    if file_key:

        buffer = "%s = %s\n" % (file_key, file_value)

        with open(rqc_file_log, "a") as myfile:
            myfile.write(buffer)

        if log: log.info("append_rqc_file: %s:%s" % (file_key, file_value))

    else:

        if log: log.warning("key or value error: %s:%s" % (file_key, file_value))

'''
Write to rqc_stats (e.g. rqc-stats.tmp) the stats_key and stats_value
@param rqc_file_log: full path to file containing key=file
@param file_key: key for the entry
@param file_value: value for the entry

'''
def append_rqc_stats(rqc_stats_log, stats_key, stats_value, log=None):

    if stats_key:

        buffer = "%s = %s\n" % (stats_key, stats_value)

        with open(rqc_stats_log, "a") as myfile:
            myfile.write(buffer)

        if log: log.info("append_rqc_stats: %s:%s" % (stats_key, stats_value))

    else:

        if log: log.warning("key or value error: %s:%s" % (stats_key, stats_value))

'''
Return the file system path to jgi-rqc-pipeline so we can use */tools and */lib

@return /path/to/jgi-rqc-pipelines
'''
def get_run_path():

    current_path = os.path.dirname(os.path.abspath(__file__))

    run_path = os.path.abspath(os.path.join(current_path, os.pardir))

    return run_path



'''
Simple function to log what happened only if it failed to run

Typical usage:
    std_out, std_err, exit_code = run_command(cmd, True, log)
    post_mortem_cmd(cmd, exit_code, std_out, std_err)

'''
def post_mortem_cmd(cmd, exit_code, std_out, std_err, log = None):

    if exit_code > 0:
        if log:
            log.error("- cmd failed: %s", cmd)
            log.error("- exit code: %s", exit_code)

        else:
            print "- cmd failed: %s" % (cmd)
            print "- exit code: %s" % (exit_code)


        if std_out:
            if log:
                log.error("- std_out: %s", std_out)
            else:
                print "- std_out: %s" % (std_out)

        if std_err:
            if log:
                log.error("- std_err: %s", std_err)
            else:
                print "- std_err: %s" % (std_err)


'''
Simple read count
- slightly different than in rqc_utility

'''

def get_read_count_fastq(fastq, log = None):

    read_cnt = 0

    if os.path.isfile(fastq):

        cat_cmd = "cat" #RQCCommands.CAT_CMD
        if fastq.endswith(".gz"):
            cat_cmd = "zcat" #RQCCommands.ZCAT_CMD
        elif fastq.endswith(".bz2"):
            cat_cmd = "bzcat" #RQCCommands.BZCAT_CMD

        cmd = "%s %s | wc -l" % (cat_cmd, fastq)
        if log:
            log.info("- cmd: %s", cmd)

        std_out, std_err, exit_code = run_command(cmd, True)




        if exit_code == 0 and std_out:
            tmp_list = std_out.split(" ")

            read_cnt = int(tmp_list[0])
            if read_cnt % 4 != 0:

                if log:
                    log.error("- bad read count for fastq, not evenly divisible by 4!  %s", read_cnt)

                read_cnt = 0

            else:
                read_cnt = read_cnt / 4
                if log:
                    log.info("- read count: %s", read_cnt)

        else:
            if log:
                post_mortem_cmd(cmd, exit_code, std_out, std_err, log)


    else:
        log.error("- fastq: %s does not exist!", fastq)

    return read_cnt



'''
Subsampling calculation
0 .. 250k reads = 100%
250k .. 25m = 100% to 1%
25m .. 600m = 1%
600m+ .. oo < 1%

July 2014 - 15 runs > 600m (HiSeq-2500 Rapid) - 4 actual libraries / 85325 seq units
- returns new subsampling rate
'''
def get_subsample_rate(read_count):

    subsample = 0
    subsample_rate = 0.01
    max_subsample = 6000000 # 4 hours of blast time

    new_subsample_rate = 250000.0/read_count
    subsample_rate = max(new_subsample_rate, subsample_rate)
    subsample_rate = min(1, subsample_rate) # if subsample_rate > 1, then set to 1

    subsample = int(read_count * subsample_rate)

    if subsample > max_subsample:
        subsample = max_subsample

    subsample_rate = subsample / float(read_count)

    return subsample_rate


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## main program



if __name__ ==" __main__":

    # unit tests

    print human_size(102192203)
    print human_size(250000000000)
    exit(0)


