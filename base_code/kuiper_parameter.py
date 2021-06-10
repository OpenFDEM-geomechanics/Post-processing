# /////////////////////////////////////////////////////////////// #
# !python3.6
# -*- coding: utf-8 -*-
# Python Script initially created on 9/7/2020
# Compiled by Aly @ Grasselli's Geomechanics Group, UofT, 2020
# Created using PyCharm // Tested on Spyder
# Current Version 06 - Dated August 21, 2018
# /////////////////////////////////////////////////////////////// #

import time
import os
import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np

# START OF EXECUTION
abs_start = time.time()

my_path = os.path.dirname(
    os.path.abspath(__file__))  # Figures out the absolute path for you in case your working directory moves around.

'''
Default MATPLOTLIB Fonts
'''

# plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams["figure.figsize"] = [5, 5]
matplotlib.rcParams['font.family'] = ['arial']
matplotlib.rcParams['font.size'] = 8

'''
TIMER FUNCTION
'''


def calc_timer_values(end_time):
    minutes, sec = divmod(end_time, 60)
    if end_time < 60:
        return "\033[1m%.2f seconds\033[0m" % end_time
    else:
        return "\033[1m%d minutes and %d seconds\033[0m." % (minutes, sec)


'''
FORMATTING OPTIONS

- TEXT COLORS
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


sens_analysis = {}

def initial_processing(frequency_scanLine):


    import matplotlib.pyplot as plt


    # A Python3 program to find if 2 given line segments intersect or not

    class Point:
        def __init__(self, x, y):
            self.x = x
            self.y = y

        # Given three colinear points p, q, r, the function checks if


    # point q lies on line segment 'pr'
    def onSegment(p, q, r):
        if ((q.x <= max(p.x, r.x)) and (q.x >= min(p.x, r.x)) and
                (q.y <= max(p.y, r.y)) and (q.y >= min(p.y, r.y))):
            return True
        return False


    def orientation(p, q, r):
        # to find the orientation of an ordered triplet (p,q,r)
        # function returns the following values:
        # 0 : Colinear points
        # 1 : Clockwise points
        # 2 : Counterclockwise

        # See https://www.geeksforgeeks.org/orientation-3-ordered-points/amp/
        # for details of below formula.

        val = (float(q.y - p.y) * (r.x - q.x)) - (float(q.x - p.x) * (r.y - q.y))
        if (val > 0):

            # Clockwise orientation
            return 1
        elif (val < 0):

            # Counterclockwise orientation
            return 2
        else:

            # Colinear orientation
            return 0


    # The main function that returns true if
    # the line segment 'p1q1' and 'p2q2' intersect.
    def doIntersect(p1, q1, p2, q2):
        # Find the 4 orientations required for
        # the general and special cases
        o1 = orientation(p1, q1, p2)
        o2 = orientation(p1, q1, q2)
        o3 = orientation(p2, q2, p1)
        o4 = orientation(p2, q2, q1)

        # General case
        if ((o1 != o2) and (o3 != o4)):
            return True

        # Special Cases

        # p1 , q1 and p2 are colinear and p2 lies on segment p1q1
        if ((o1 == 0) and onSegment(p1, p2, q1)):
            return True

        # p1 , q1 and q2 are colinear and q2 lies on segment p1q1
        if ((o2 == 0) and onSegment(p1, q2, q1)):
            return True

        # p2 , q2 and p1 are colinear and p1 lies on segment p2q2
        if ((o3 == 0) and onSegment(p2, p1, q2)):
            return True

        # p2 , q2 and q1 are colinear and q1 lies on segment p2q2
        if ((o4 == 0) and onSegment(p2, q1, q2)):
            return True

        # If none of the cases
        return False

    # This code is contributed by Ansh Riyal



    id = {}

    import csv

    # Sample Dimensions
    w = 54
    h = 108
    # File Name
    # Load File and get X/Y as x0,x1,y0,y1
    f_name = '/home/aly/Desktop/sourcejoints.csv'
    with open(f_name)as csv_file:
        reader = csv.reader(csv_file,  skipinitialspace=True, )
        next(reader)
        for line in reader:
            id[line[1]] = [float(line[2]), float(line[5]), float(line[3]), float(line[6])]
    # print(id)

    # First figure that shows the tracelines and fractures
    fig, ax = plt.subplots()
    ax.set_xlabel("Sample Width (mm)")
    ax.set_ylabel("Sample Heigth (mm)")
    ax.set_title("Frequency of Scan Line [%0.2f]" % frequency_scanLine)
    # Plot fractures
    for k,v in id.items():
        ax.plot([v[0], v[1]], [v[2], v[3]])

    l_kuiper = []
    l_x, l_y = [], []
    counter = 0

    # Plot trace lines on fractures
    for i in np.arange(0, h, frequency_scanLine):
        ax.plot([0, w], [i, i])
        p1 = Point(0, i)
        q1 = Point(w, i)
        for k, v in id.items():
            # l1 = LineString([(v[0], v[2]), (v[1], v[3])])
            p2 = Point(v[0], v[2])
            q2 = Point(v[1], v[3])

            if doIntersect(p1, q1, p2, q2):
                # Cumulative Counter
                counter += 1

        l_kuiper.append([i, counter])
        l_x.append(i)
        l_y.append(counter)


    # print("L Kuiper", l_kuiper)
    fig1, ax1 = plt.subplots()
    # Plot the Stepping cumulative curve
    ax1.step(l_x, l_y, linestyle=":", label='Cumulative graph', alpha=0.75, color='blue')
    # ax1.plot(l_x, l_y, '--', color='grey', alpha=0.3)

    myradians=[]
    d_plus = 0
    d_minus= 0

    # Draw UNIFORM line
    ax1.plot([0, max(l_x)], [0, max(l_y)], linestyle="--", label='Uniform', color='orange')
    # Calculate angle of UNIFORM line
    slope = max(l_y) / max(l_x)
    uni_ang = math.degrees(math.atan2(max(l_y),  max(l_x)))
    print(bold_text("\nFrequency of Scan Line %.2f" % frequency_scanLine))
    print("Angle of Uniform line %.2f degrees" % uni_ang)
    for i in l_kuiper:
        ax1.scatter(i[0], i[1], marker='x', color='black')
        # Get co-ord of uniform line WRT to line
        # Calculate the dy between the two
        x1= i[0]
        y1= i[0] * math.tan(math.radians(uni_ang))
        dy = y1 - i[1]
        # print("Trace %02d had P_{10} of %02d with distance to UNIFORM line %.2f" % (i[0], i[1], dy))
        ang = math.atan2(i[1] - 0, i[0] - 0)
        myradians.append(ang)
        if dy > d_plus:
            d_plus = dy
            idx_max = i
        if dy < d_minus:
            d_minus = dy
            idx_min = i


    def plot_point(point, angle, length):
        '''
        point - Tuple (x, y)
        angle - Angle you want your end point at in degrees.
        length - Length of the line you want to plot.

        Will plot the line on a 10 x 10 plot.
        '''

        # unpack the first point
        x, y = point

        # find the end point
        endy = length * math.sin(math.radians(angle)) + y
        endx = length * math.cos(math.radians(angle)) + x
        return [x, endx], [y, endy]

    def draw_parallael_line(coord, slope, lab, col):
        m, p = coord
        axes = plt.gca()
        x_val = np.array(axes.get_xlim())
        y_val = np.array(slope * (x_val - m) + p)
        plt.plot(x_val, y_val, color=col, linestyle="--", label=lab)

    print("D Values:\n\tMax: %.2f\n\tMin: %.2f\n\tCum Total: %d " % (d_plus, d_minus, counter))
    ax1.set_xlim(0, h)
    ax1.set_ylim(0, counter)
    ax1.set_xlabel("Distance (mm)")
    ax1.set_ylabel("Cumulative Values (-)")


    draw_parallael_line(idx_max, slope, "D-", 'red')
    draw_parallael_line(idx_min, slope, 'D+', 'green')

    ax1.legend()

    s_hetero = (abs(d_plus) + abs(d_minus)) / counter
    print("Hetro Ratio %0.5f" % s_hetero)
    ax1.set_title("Kuiper Test %0.5f [%.2f]" % (s_hetero, frequency_scanLine))
    print("The result is a number between 0, perfect strain homogeneity, and 1, maximum strain heterogeneity, which is possible if all strain is manifested in a single fracture ")

    from scipy import stats
    print(stats.kstest(l_y, 'norm'))
    sens_analysis[frequency_scanLine] = s_hetero
    plt.tight_layout()
    # fig1.tighlayout()
    print(os.path.join(my_path, "Sample_ScanLine_%s.png" % str(frequency_scanLine).replace(".","_")))
    fig.savefig(os.path.join(my_path, "Sample_ScanLine_%s.png" % str(frequency_scanLine).replace(".","_")), dpi=600)
    fig1.savefig(os.path.join(my_path, "Kuiper_Test_%s.png" % str(frequency_scanLine).replace(".","_")), dpi=600)
    plt.show()
    return sens_analysis

def plot_bar_chart_from_Dict(a_dictionary):
    keys = a_dictionary.keys()
    values = a_dictionary.values()

    plt.bar(range(len(a_dictionary)), list(a_dictionary.values()), align='center')
    plt.xticks(range(len(a_dictionary)), list(a_dictionary.keys()))
    plt.title('Sensitivity Analysis')
    plt.xlabel("Distance of Scan Line (mm)")
    plt.ylabel("Kuiper Test Value")
    plt.savefig(os.path.join(my_path, "Kuiper_Test_Sensitivity_Analysis.png"), dpi=600)
    plt.show()

'''
MAIN MODULE

- Returns total time and Error on user termination.
'''

if __name__ == "__main__":
    try:
        for freq in [0.25, 0.5, 1, 2, 5, 10, 15]:
            initial_processing(freq)
        plot_bar_chart_from_Dict(initial_processing(freq))
        print("\nTotal Execution time: \033[1m%s\033[0m\n" % calc_timer_values(time.time() - abs_start))
    except KeyboardInterrupt:
        exit("TERMINATED BY USER")
