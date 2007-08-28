#!/usr/bin/python
# -*- coding: UTF-8 -*-
# 
# Calculating genotype frequencies and storing them in the database.
#
#
# Copyright (C) 2007 Maciej Bliziński
# 
# This script is free software distributed WITHOUT ANY WARRANTY. 
# You can redistribute it and/or modify it under the terms of the
# GNU General Public License, 
# version 2 or later, as published by the Free Software Foundation. 
# See the file COPYING for details.
# 
 
import psycopg
import re
import itertools
import time
import sys
from userdata import group_list

class GenotypeCounter(object):
    batch_size = 1000
    def __init__(self, conn):
        self.conn = conn
        self.inscur = conn.cursor()
        self.row_updated = False
    def go(self):
        curs = self.conn.cursor()
        # curs.execute("""SELECT count(locus_id || data_set_name) FROM locus_population;""")
        # total_rows = curs.fetchone()[0]

        # Hard-coding this for speed:
        curs.execute("SELECT max(id) FROM locus_population;")
        total_rows = curs.fetchone()[0]
        self.batch_size = total_rows / 100 + 1
        print "Total rows:", total_rows
        # Batches
        # limit is stable, say, 50
        # offset is 0, 50, 100, and so on.
        batches = xrange(0, total_rows, self.batch_size)
        qry = """
SELECT data_set_name, locus_id, indivs, g1freq
FROM locus_population
-- WHERE g1freq IS NOT NULL
-- ORDER BY data_set_name ASC, locus_id ASC
WHERE id >= %s AND id < (%s + %s)
AND g1freq IS NULL
;"""
        for offset in batches:
            self.row_updated = False
            sys.stdout.write("querying with offset: %s ... " % offset)
            sys.stdout.flush()
            curs.execute(qry, (offset, offset, self.batch_size))
            sys.stdout.write("processing rows... ")
            sys.stdout.flush()
            while True:
                row = curs.fetchone()
                if not row:
                    break
                self.__cnt_gt(row[0], row[1], row[2])

            if self.row_updated == True:
                sys.stdout.write("committing... ")
                sys.stdout.flush()
                self.conn.commit()
                sys.stdout.write("sleeping... ")
                sys.stdout.flush()
                time.sleep(5)
            sys.stdout.write("done.\n")
            sys.stdout.flush()
    def __cnt_gt(self, data_set_name, locus_id, indivs):
        gts = re.findall(r'[0-9\.]+', indivs.strip())
        grouped = group_list(gts, 2)
        a1 = gts.count('1')
        a2 = gts.count('2')
        gametes = len(gts)
        assert (a1 + a2) == gametes, "Problem with a1+a2 vs gametes"
        # print grouped
        # Count
        d = {
                ('1', '1'): 0,
                ('1', '2'): 0,
                ('2', '1'): 0,
                ('2', '2'): 0, }
        for elem in grouped:
            d[elem] += 1
        # map 2,1 to 1,2, since we don't want to distinguish
        d[('1', '2')] += d[('2', '1')]
        d.pop(('2', '1'))
        self.inscur.execute("""
UPDATE locus_population
SET g1freq = %s, g2freq = %s, g3freq = %s,
gametes = %s, a1 = %s, a2 = %s
WHERE data_set_name = %s AND locus_id = %s;""",
        (d[('1', '1')], d[('1', '2')], d[('2', '2')],
            gametes, a1, a2,
            data_set_name, locus_id,))
        self.row_updated = True

if __name__ == '__main__':
	conn = psycopg.connect("dbname='genepi'")
	gc = GenotypeCounter(conn)
	gc.go()

