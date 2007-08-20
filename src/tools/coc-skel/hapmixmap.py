#!/usr/bin/python
# -*- coding: UTF-8 -*-
# 
# Classes for HTML generator for hapmixmap results.
#
# This file is part of Genepi, genetic data analysis software.
# Copyright (C) 2007 Maciej Bliziński
# 
# Genepi is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# Genepi is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Genepi; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

import os
import os.path
import re
import shutil
from time import localtime, strftime
import colorsys

def slugify(s):
    return '-'.join(re.findall(r'\w+', s.lower()))

class TextFile(object):
    def __init__(self, path, fn, out_dir):
        self.path = path
        self.fn = fn
        self.out_dir = out_dir
    def get_dst_bn(self):
        return os.path.basename(self.get_dst_fn())
    def get_dst_fn(self):
        return os.path.join(
                self.out_dir,
                slugify(self.path) + "-" + self.fn)
    def get_src_fn(self):
        return os.path.join(self.path, self.fn)
    def get_html(self):
        return "<a href=\"%s\">%s</a>" % (
                self.get_dst_bn(),
                self.fn)

class PostScriptImage(object):
    def __init__(self, fn, out_dir, group_id):
        assert fn[-3:] == '.ps'
        self.base = fn[:-3]
        self.out_dir = out_dir
        self.group_id = group_id
        assert os.path.exists(fn), "Image file missing."
        # TODO: Check the timestamps to find out if the file needs to be
        # updated.
        if not os.path.exists(self.get_gif_fn()):
            command = self.get_convert_command(self.get_gif_fn(), False)
            print command; os.popen(command).readlines()
        if not os.path.exists(self.get_gif_thumbnail_fn()):
            command = self.get_convert_command(self.get_gif_thumbnail_fn(), True)
            print command; os.popen(command).readlines()
        shutil.copy(self.get_gif_fn(), self.get_dst_gif_fn())
        shutil.copy(self.get_gif_thumbnail_fn(), self.get_dst_gif_thumbnail_fn())
    def get_convert_command(self, dst_fn, resize = False):
        resize_par = ""
        if resize:
            resize_par = "-resize 180x120"
        return "convert -rotate 90 -delay 200 %s %s %s" % \
                (self.get_ps_fn(), resize_par, dst_fn)
    def get_gif_thumbnail_fn(self):
        return "%s-thumb.gif" % self.base
    def get_gif_fn(self):
        return "%s.gif" % self.base
    def get_ps_fn(self):
        return "%s.ps" % self.base
    def get_dst_gif_fn(self):
        return os.path.join(self.out_dir, "%s-%s" % \
                (self.get_identity_prefix(), os.path.basename(self.get_gif_fn())))
    def get_dst_gif_thumbnail_fn(self):
        return os.path.join(self.out_dir, "%s-%s" % \
                (self.get_identity_prefix(), os.path.basename(self.get_gif_thumbnail_fn())))
    def get_identity_prefix(self):
        return slugify(self.group_id)
    def get_html(self):
        return """
        <a href="%s" title="%s"><img class="chart" src="%s" alt="%s"/></a>
        """ % (os.path.basename(self.get_dst_gif_fn()),
                            "Zoom in " + os.path.basename(self.base),
                            os.path.basename(self.get_dst_gif_thumbnail_fn()),
                            os.path.basename(self.base),)

class LinearScaler(object):
    def __init__(self, src_low, src_up, targ_low, targ_up):
        self.src_low  = float(src_low)
        self.src_up   = float(src_up)
        self.targ_low = float(targ_low)
        self.targ_up  = float(targ_up)
    def scale(self, x):
        ans = 0.0
        if x <= self.src_low:
            ans = self.targ_low
        elif x >= self.src_up:
            ans = self.targ_up
        else:
            ans = (x - self.src_low) / (self.src_up - self.src_low) \
                * (self.targ_up - self.targ_low) + self.targ_low
        return ans

class ColorCalculator(object):
    def __init__(self):
        self.scalers = []
        # Hue
        self.scalers.append(LinearScaler(0.7, 0.75, 0.3, 0.7))
        # Saturation
        self.scalers.append(LinearScaler(0.7, 0.75, 0.0, 0.3))
        # Value
        self.scalers.append(LinearScaler(0.7, 0.75, 0.97, 0.93))
        self.norm_to_byte_scaler = LinearScaler(0, 1, 0, 255)
    def get_rgb(self, value):
        h, s, v = [x.scale(value) for x in self.scalers]
        val_list = colorsys.hsv_to_rgb(h, s, v)
        val_list = [self.norm_to_byte(x) for x in val_list]
        return ", ".join([str(int(x)) for x in val_list])
    def norm_to_byte(self, value):
        return int(self.norm_to_byte_scaler.scale(value))

def cartesian_product(lists, previous_elements = []):
    """Generates a cartesian product of given list of lists."""
    assert type(lists) == list
    assert len(lists) >= 1
    assert type(lists[0]) == list
    if len(lists) == 1:
        for elem in lists[0]:
            yield previous_elements + [elem, ]
    else:
        for elem in lists[0]:
            for x in cartesian_product(lists[1:], previous_elements + [elem, ]):
                yield x

class Parameter(object):
    def __init__(self, long_name, short_name, param_name, values):
        assert type(values) == list
        self.short_name = short_name
        self.long_name = long_name
        self.param_name = param_name
        self.values = values
    def __str__(self):
        return self.long_name
    def __repr__(self):
        return "<P: %s>" % (self.long_name)

class ParameterSet(object):
    """Defines a set of parameters for an analysis.
    Program will be run for every element of the cartesian
    product of the parameter value sets."""
    def __init__(self):
        self.parameters = []
    def add_parameter(self, parameter):
        assert parameter.__class__.__name__ == 'Parameter'
        self.parameters.append(parameter)
    def get_all_combinations(self):
        """Return a list of lists of pairs, where first element is
        a parameter and the second is parameter's value. """
        # return [(self.parameters, x) for x in cartesian_product(
        # [x.values for x in self.parameters])]
        cpp_pairs = []
        for cpe in cartesian_product([x.values for x in self.parameters]):
            assert len(cpe) == len(self.parameters)
            counter = 0L
            pairs = []
            for idx in cpe:
                pairs.append((self.parameters[counter], cpe[counter]))
                counter += 1
            cpp_pairs.append(pairs)
        return cpp_pairs

class HtmlAwareError(object):
    def __float__(self):
        return 0.0
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return '<acronym title="%s">Error</a>' % str(self.msg)

class TextTable(object):
    def __init__(self):
        self.headers = []
        self.tbody = []
        self.style = {}
    def read_file(self, file_name):
        self.file_name = file_name
        col_splitter = re.compile(r'\s+')
        tcl = open(self.file_name).readlines()
        self.headers = col_splitter.split(tcl[0])
        for row in [col_splitter.split(x) for x in tcl[1:]]:
            self.add_row(row)
    def get_html(self):
        ret = "<table>\n"
        ret += "<tr>\n"
        ret += " ".join([("<th>%s</th>\n" % x) for x in self.headers])
        ret += "</tr>\n"
        for row in self.tbody:
            ret += "<tr>\n"
            for elem in row:
                assert type(elem) == dict
                if elem.has_key('style'):
                    style = " style=\"%s\"" % elem['style']
                else:
                    style = ""
                ret += ' <td class="numeric"%s>%s</td>\n' % (style, elem['content'])
            ret += "</tr>\n"
        ret += "</table>\n"
        return ret
    def set_headers(self, headers):
        assert type(headers) == list
        self.headers = headers
    def add_row(self, row):
        assert type(row) == list
        """If it's the old-style value elements, convert them
        to dict-type elements (with style)."""
        new_row = []
        for elem in row:
            if type(elem) != dict:
                tmp = {}
                tmp['content'] = elem
                elem = tmp
            new_row.append(elem)
        self.tbody.append(new_row)

class HapmixmapLog(object):
    def __init__(self, file_name):
        self.file_name = file_name
        try:
            self.content_lines = open(self.file_name, "r").readlines()
        except IOError, e:
            print "Could not open", self.file_name
            self.content_lines = None
    def get_execution_time(self):
        try:
            time_re = re.compile("Program finished")
            whitespace_re = re.compile("\s+")
            slash_re = re.compile("/")
            for line in self.content_lines:
                if time_re.search(line):
                    token_list = whitespace_re.split(line)
                    date_token = token_list[4]
                    time_token = token_list[3]
                    date_list = slash_re.split(date_token)
                    return "%04d-%02d-%02d %s" % (
                            int(date_list[2]),
                            int(date_list[1]),
                            int(date_list[0]),
                            token_list[3])
        except:
            return "date unknown"

class PopulationResults:
    """Results for a single population."""
    def __init__(self, population, base_dir, slug):
        self.population = population
        self.base_dir = base_dir
        self.slug = slug
        self.images = []
        self.parameters = {}
        self.log = None
    def get_log(self):
        if not self.log:
            self.log = HapmixmapLog(self.get_log_file_name())
        return self.log
    def get_html(self):
        energy_plot_template = open("template/energy-plot.html", "r").read()
        energy_plots = ''
        # thumb_re = re.compile(r'.*200.*thumb.gif')
        thumb_re = re.compile(r'.*thumb.gif')
        for image in self.images:
            if not thumb_re.search(image):
                continue
            energy_plots += energy_plot_template % ({
                'population': self.population,
                'slug': self.slug,
                'file_name': image,
                'big_file_name': image.replace('-thumb.gif', '.gif'), })
        return open("template/population-template.html").read() % {
            'population': self.population,
            'posterior_quantiles': self.get_posterior_quantiles(), 
            'slug': self.slug,
            'energy_plots': energy_plots, }
    def get_results_dir(self):
        return os.path.join(self.base_dir, "hapmap", self.population, "Chr22Results6States2")
    def get_log_file_name(self):
        return os.path.join(self.get_results_dir(), "logfile.txt")
    def get_coefficient_of_constraint(self):
        try:
            return "%05.03f" % float(open(os.path.join(
                self.get_results_dir(),
                "mean-coefficient-of-constraint.txt")).read())
        except IOError, e:
            return HtmlAwareError(e)
    def get_posterior_quantiles(self):
        t = TextTable()
        try:
            t.read_file(os.path.join(self.get_results_dir(), "PosteriorQuantiles.txt"))
            return t.get_html()
        except IOError, e:
            return HtmlAwareError(e)
    def get_files(self, names):
        assert type(names) == list
        for root, dirs, files in os.walk(self.get_results_dir()):
            for file in files:
                print file
    def get_execution_time(self):
        return self.get_log().get_execution_time()

class ThreePopulationResults:
    populations = ['Afr', 'Asian', 'Eur']
    def __init__(self, slug, base_dir):
        self.slug = slug
        self.base_dir = base_dir
        self.populations_results = []
        for population in self.populations:
            self.populations_results.append(
                    PopulationResults(
                        population,
                        os.path.join(self.base_dir, self.slug),
                        self.slug))
    def get_coefficient_of_constraint_for_3_populations(self):
        t = TextTable()
        t.set_headers(self.populations)
        the_only_row = []
        color_c = ColorCalculator()
        for population in self.populations_results:
            coc = population.get_coefficient_of_constraint()
            coc_dict = {}
            coc_dict['content'] = coc
            try:
                coc_dict['style'] = "background-color: rgb(%s);" % \
                        (color_c.get_rgb(float(coc)))
            except AttributeError, e:
                print "AttributeError:", e
            the_only_row.append(coc_dict)
        t.add_row(the_only_row)
        return t.get_html()
    
    def get_text_coc_3_pop(self):
        ret = ""
        coc_list = []
        for population in self.populations_results:
            coc_list.append(population.get_coefficient_of_constraint())
        return ", ".join([str(x) for x in coc_list])

    def get_html(self):
        populations_html = reduce(lambda x, y: x + y, map(lambda x: x.get_html(), self.populations_results), '')
        return open("template/three-populations-template.html", "r").read() % {
            'coefficient_of_constraint': self.get_coefficient_of_constraint_for_3_populations(),
            'populations': populations_html, 
            'slug': self.slug, }

    def copy_file(self, dest_dir, file_name):
        for population in self.populations_results:
            # copy the argsfile
            try:
                shutil.copyfile(
                    os.path.join(population.get_results_dir(), file_name),
                    os.path.join(dest_dir, "%s-%s-%s" % (
                        self.slug, population.population, file_name)))
            except IOError, e:
                print e
    def copy_image_by_regex(self, dest_dir, file_search):
        for population in self.populations_results:
            for root, dirs, files in os.walk(population.get_results_dir()):
                for file in files:
                    if file_search.match(file):
                        new_name = "%s-%s-%s" % (
                                self.slug, population.population, file)
                        shutil.copyfile(
                                os.path.join(root, file),
                                os.path.join(dest_dir, new_name))
                        population.images.append(new_name)
    def get_execution_time(self):
        return self.populations_results[0].get_execution_time()

class AnalysisResultsGroup:
    """Class representing results from a set of analysis, as defined
    by a `analysis.txt' file."""
    def __init__(self, main_dir, dest_dir):
        assert type(main_dir) == str
        assert type(dest_dir) == str
        self.main_dir = main_dir
        self.analyses = {}
        self.dest_dir = dest_dir
        self.html_template = open("template/results-template.html").read()
    def add_src_directory(self, src_dir):
        assert type(src_dir) == str
        self.analyses[src_dir] = ThreePopulationResults(src_dir, self.main_dir)
    def write_html(self):
        # print "AnalysisResultsGroup::write_html()"
        exec_times, by_exec_time = self.get_results_by_exec_time()
        exec_times_len = len(exec_times)
        for exec_time_no in range(exec_times_len):
            exec_time = exec_times[exec_time_no]
            tpl = by_exec_time[exec_time]
            src_dir = tpl.slug
            prev_slug = by_exec_time[exec_times[(exec_time_no - 1) % exec_times_len]].slug
            next_slug = by_exec_time[exec_times[(exec_time_no + 1) % exec_times_len]].slug

            # copy_file_by_regex collects information about images to be
            # used during HTML generation.
            tpl.copy_image_by_regex(
                    self.dest_dir,
                    # re.compile(r"EnergyTracePlot-[0-9].*\.gif"))
                    re.compile(r"EnergyTracePlot.*\.gif"))

            open(os.path.join(self.dest_dir, src_dir + ".html"), "w").write(
                    self.html_template % {
                        'html_body': tpl.get_html(),
                        'html_title': src_dir,
                        'date': tpl.get_execution_time(),
                        'prev_slug': prev_slug,
                        'next_slug': next_slug, })
            tpl.copy_file(self.dest_dir, "args.txt")
            tpl.copy_file(self.dest_dir, "logfile.txt")

        # Write index HTML file
        fh = open(os.path.join(self.dest_dir, "index.html"), "w")
        fh.write(open(os.path.join("template", "index.html")).read() % {
            'date': '', # FIXME
            'toc': '', # FIXME
            'index_table': self.get_index_table(), })
        fh.close()
        shutil.copyfile("template/hapmixmap.css",
                os.path.join(self.dest_dir, "hapmixmap.css"))
    def get_results_by_exec_time(self):
        by_exec_time = {}
        for key in self.analyses.keys():
            tpl = self.analyses[key]
            by_exec_time[tpl.get_execution_time()] = tpl
        exec_times = by_exec_time.keys()
        exec_times.sort()
        return exec_times, by_exec_time
    def get_index_table(self):
        ret = '<table>\n'
        exec_times, by_exec_time = self.get_results_by_exec_time()
        for exec_time in exec_times:
            tpl = by_exec_time[exec_time]
            ret += "<tr>\n"
            ret += "<td><a href=\"%s.html\">%s</a></td>\n" % (tpl.slug, tpl.slug)
            ret += "<td>%s</td>\n" % (tpl.get_execution_time())
            ret += "<td>%s</td>\n" % (tpl.get_coefficient_of_constraint_for_3_populations())
            ret += "</tr>"
        ret += '<table>\n'
        return ret
    def get_index_csv(self):
        ret = ""
        exec_times, by_exec_time = self.get_results_by_exec_time()
        for exec_time in exec_times:
            tpl = by_exec_time[exec_time]
            ret += "%s, %s\n" % (tpl.slug, tpl.get_text_coc_3_pop())
        return ret
    def write_csv(self):
        fh = open(os.path.join(self.dest_dir, "results.csv"), "w")
        fh.write(self.get_index_csv())
        fh.close()

