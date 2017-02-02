# =====================================================================================
# | A new algorithm for spatial clustering of the line segments
# | as a tool for lineament extraction.
# |
# | Ondrej Kaas, 2016
# | Faculty of Applied Sciences,University of West Bohemia, Pilsen, Czech Republic
# =====================================================================================

# GLOBAL MODULES =============================
import math
import time
import random
import argparse

# ============================================

# OWN STRUCTURES =============================
class Vertex:
    def __init__(self, coords, i):
        self.coords = []
        for c in coords[:-1]:
            self.coords.append(int(c))

        self.clients = []
        self.to_fac = 0.0
        self.acumulator = 0.0
        self.index = i

    # def info(self):
    #     if len(self.clients) > 0:
    #         print '{}'.format(self.clients)
    #     print self.coords, sep=', '

    def remove_client(self, index):
        for i in range(0, len(self.clients)):
            if self.clients[i].index == index:
                self.clients.remove(self.clients[i])
                break


# ============================================

sX = 1
sY = 2
eX = 3
eY = 4
aZ = 5
LEN = 6

WEIGHTS = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
FACILITY_COST = 1

BORDER_X = -1
BORDER_Y = -1
BORDER_AZ = -1
OUTPUT_FILE = ""

VERTICES = []
VERTICES_LEN = -1
VERTICE_FOR_PERFORM = []
FACILITIES = []

REASSIGN_VERTICES = []
FACILITIES_CLOSE = []

INFINITE = float('inf')
INFINITE_NEG = -float('inf')

MIN_DIST = INFINITE



def read_data(input_file):
    """
    Read data from input file and save it into list of vertices
    :param input_file: path to input file
    :return: list of vertices
    """
    f = open(input_file, 'r')

    # skip first line
    f.readline()

    index = -1

    for line in f:
        index += 1
        VERTICES.append(Vertex(line.split(';'), index))

    return VERTICES


def weighted_distance(a, b):
    """
    Calculate weighted distance between two vertices
    :param a: coord of vertex A
    :param b: coord vertex B
    :return: distance
    """
    diff_sqr_sum = 0
    dim = 4

    for i in range(0, dim):
        d = a[i] - b[i]
        diff_sqr_sum += d * d

    return math.sqrt(diff_sqr_sum)


def reassign(vertex_index):
    """
    Perform reassign of vertex at index "vertex_index" for better cluster solution
    :param vertex_index: index of vertex its connection will be changed
    :return: None
    """
    # facility where to reassign
    facil_to_assign = -1

    # create the new facility
    if VERTICES[vertex_index].clients is None:
        # create a new facility
        facil_to_assign = VERTICES[vertex_index]
        FACILITIES.insert(vertex_index, facil_to_assign)

    # get the new facility
    new_facility = VERTICES[vertex_index]

    # perform re-assignments
    for index in range(0, len(VERTICE_FOR_PERFORM)):
        # remove from current facility
        facilRemoveFrom = VERTICES[index]
        facilRemoveFrom.remove_client(index)

        # assign to new facility
        facil_to_assign.clients.append(index)

    # perform facility closures
    for index_fac in range(0, len(FACILITIES_CLOSE)):
        # reassign all it's points to the newly created facility (may be precomputed)
        # facility vertex itself will be re-assigned too
        for index_clients in range(0, VERTICES[FACILITIES_CLOSE[index_fac].index].clients):
            # removing vertex from current facility not neccessary
            # we will close the facility anyway
            # removing even not possible because we iterate over the list

            # set new facility
            dist = VERTICES[index].WeightedDistance(new_facility)
            VERTICES[index].to_fac = dist

            # assign to new facility
            facil_to_assign.clients.append(index)

        # remove facility from list
        FACILITIES.remove(index_fac)

        # vertex is no longer a facility
        VERTICES[FACILITIES_CLOSE[index_fac].index].clients = None


def gain(vertex_index):
    """
    Function which determine income in overall clustering cost if vertex at index "vertex_index" change its connection
    :param vertex_index: index of vertex its new connection will be investigated
    :return: potential improvement by changing current connection of vertex at index "vertex_index"
    """
    gain = 0

    # vertex is facility > pay facility cost
    if VERTICES[vertex_index].clients is not None:
        gain -= FACILITY_COST

    for i in range(0, VERTICES_LEN):
        dist_gain = VERTICES[vertex_index].to_fac - distance(VERTICES[i], VERTICES[vertex_index])

        if dist_gain > 0:
            # reassign_vertices.append(i)
            gain += dist_gain
        else:
            VERTICES[vertex_index].acumulator += dist_gain

    for i in range(0, len(FACILITIES)):
        if FACILITIES[i].acumulator + FACILITY_COST and i != vertex_index:
            FACILITIES_CLOSE.append(i)

            gain += FACILITIES[i].acumulator + FACILITY_COST

        FACILITIES[i].acumulator = 0.0
        FACILITIES[i].acumulator = 0.0

    return gain


def print_result(file_name):
    """
    Write clustering input into file
    :param file_name: path to output file
    :return: None
    """
    write = open(file_name, 'w')

    for v in VERTICES:
        write.write("> ")
        for coords in v.coords:
            write.write("%d;" % coords)

        write.write("\n")

        if len(v.clients) > 0:
            for cl in v.clients:
                write.write("* ")
                for cl_coord in cl.coords:
                    write.write("%d;" % cl_coord)
                write.write("\n")

        write.write("\n")


def distance(a, b):
    """
    Calculate distance between two lines A and B
     > if A is in buffer zone of B => return weighted distance between them
     > else return infinity
    :param a:
    :param b:
    :return:
    """
    if in_buffer(a.coords, b.coords):
        return 0.0
    else:
        return weighted_distance(a.coords, b.coords)


def in_buffer(a, b):
    """
    Check if line B is in buffer zone of A
    :param a: coord of line A
    :param b: coord of line B
    :return: False - line B is NOT in buffer zone of A
             True  - line B is in buffer zone of A
    """
    # directional vector for A
    dir_a = [a[eX] - a[sX], a[eY] - a[sY]]

    # lenght of A
    len_a = math.sqrt(dir_a[0] * dir_a[0] + dir_a[1] * dir_a[1])

    # coeficient for shift in X
    shift_x = BORDER_X / len_a

    # coeficient for shift in Y
    shift_y = BORDER_Y / len_a

    # checking azimuth
    if b[aZ] < BORDER_AZ or b[aZ] > (180 - BORDER_AZ):
        blue_min = b[aZ] - BORDER_AZ

        if blue_min < 0:
            blue_min += 180

        blue_max = b[aZ] + BORDER_AZ

        if blue_max > 180:
            blue_max -= 180

        # checking if azimuth is in tolerance
        if not (((a[aZ]) >= 0 and a[aZ] < blue_max) or (blue_min < a[aZ] < 180)):
            return False
    else:
        if a[aZ] < b[aZ]:
            if (b[aZ] - a[aZ]) > BORDER_AZ:
                return False
        else:
            if (a[aZ] - b[aZ]) > BORDER_AZ:
                return False

    # Name of particular tested lines around start line (A)
    #
    #              3
    #      ----------------
    #     |        |       |
    # 4   | s --------- e  | 2
    #     |        |       |
    #      ----------------
    #              1
    #

    # shifting border line "behind" the starting point
    four = [a[sX] + (-shift_x) * dir_a[0], a[sY] + (-shift_x) * dir_a[1]]

    # normal vector of 4th line
    c_four = -(four[0] * dir_a[0]) - four[1] * dir_a[1]

    # general equation is > ax + by + c = 0
    # based on mark of value c after appointment for 4th line is possible decide where teste line is
    # positive mark means inside, negative outside
    if (b[sX] * dir_a[0] + b[sY] * dir_a[1] + c_four) < 0 or (b[eX] * dir_a[0] + b[eY] * dir_a[1] + c_four) < 0:
        return False

    # normal vector of 2th line
    second = [a[eX] + shift_x * dir_a[0], a[eY] + shift_x * dir_a[1]]

    # value of c for 2th line
    c_second = -(second[0] * dir_a[0]) - second[1] * dir_a[1]

    # positive mark means outside, negative inside
    if (b[sX] * dir_a[0] + b[sY] * dir_a[1] + c_second) > 0 or (b[eX] * dir_a[0] + b[eY] * dir_a[1] + c_second) > 0:
        return False

    # shifting border line "forward" the starting point
    first = [a[sX] + (-shift_y) * (-dir_a[1]), a[sY] + (-shift_y) * dir_a[0]]

    # value of c for first line
    c_first = -(first[0] * -dir_a[1]) - (first[1] * dir_a[0])

    # positive mark means inside, negative outside
    if (b[sX] * -dir_a[1] + b[sY] * dir_a[0] + c_first) < 0 or (b[eX] * -dir_a[1] + b[eY] * dir_a[0] + c_first) < 0:
        return False

    # shifting border 3th line
    third = [a[sX] + shift_y * (-dir_a[1]), a[sY] + shift_y * (dir_a[0])]

    # value of c for 3th line
    c_third = -(third[0] * (-dir_a[1])) - third[1] * dir_a[0]

    # positive mark means outside, negative means inside
    if (b[sX] * -dir_a[1] + b[sY] * dir_a[0] + c_third) > 0 or (b[eX] * -dir_a[1] + b[eY] * dir_a[0] + c_third) > 0:
        return False

    # if all previous conditions does not be hit we can be sure line A and B are in buffer zone
    return True


def compute_bounding_box():
    """
    Calculate minimum and maximal value in particular axis
    :return: None
    """
    global border_max
    global border_min

    border_max = [-INFINITE_NEG] * len(VERTICES[0].coords)
    border_min = [INFINITE] * len(VERTICES[0].coords)

    # sum up the squares of the box sizes in all dimensions
    for i in range(0, len(VERTICES)):
        # multiply the size by the weight

        for dim in range(0, len(VERTICES[0].coords)):
            if border_max[dim] < VERTICES[i].coords[dim]:
                border_max[dim] = VERTICES[i].coords[dim]

            if border_min[dim] > VERTICES[i].coords[dim]:
                border_min[dim] = VERTICES[i].coords[dim]


def closest_facility(vertex):
    """
    Find the closest facility to given vertex
    :param vertex: vertex index
    :return: None - if no facility is not in buffer zone of given vertex
             Facility - the closest facility
    """
    global MIN_DIST

    MIN_DIST = Ellipsis
    ret = -1

    for facility in range(0, int(len(FACILITIES))):

        in_buff = in_buffer(VERTICES[facility].coords, vertex.coords)

        if in_buff:
            return FACILITIES[ret]

    return None


def generate_initial_solution():
    """
    Generate coarse solution of clustering
    :return:
    """
    global MIN_DIST
    global border_min
    global border_max

    random.shuffle(VERTICES)

    # first vertex is always facility
    FACILITIES.append(VERTICES[0])

    for i in range(1, VERTICES_LEN):
        # start_closest = time.clock()

        fac = closest_facility(VERTICES[i])

        if fac is None:
            VERTICES[i].clients.append(VERTICES[i])
            FACILITIES.insert(i, VERTICES[i])
        else:
            fac.clients.append(VERTICES[i])


def perform_iterations(iter):
    """
    Try to perform improvement which helps in overall clustering cost
    :param iter: number of improvement iteration
    :return:
    """
    for i in range(0, iter):
        vertex_index = random.randint(0, VERTICES_LEN - 1)

        gain = 0

        # vertex is facility > pay facility cost
        if VERTICES[vertex_index].clients is not None:
            gain -= FACILITY_COST

        for i in range(0, VERTICES_LEN):
            dist_gain = VERTICES[vertex_index].to_fac - distance(VERTICES[i], VERTICES[vertex_index])

            if dist_gain > 0:
                gain += dist_gain
            else:
                VERTICES[vertex_index].acumulator += dist_gain

        for i in range(0, len(FACILITIES)):
            if FACILITIES[i].acumulator + FACILITY_COST and i != vertex_index:
                gain += FACILITIES[i].acumulator + FACILITY_COST

                if gain > 0:
                    reassign(vertex_index)


def clustering():
    """
    Main function of facility location clustering
    :return:
    """
    start_total = time.clock()

    # find bounding box
    compute_bounding_box()

    generate_initial_solution()
    perform_iterations(int(math.log10(VERTICES_LEN / 10)))

    print "Hello world"
    print "Lines:\t\t\t {}".format(VERTICES_LEN)
    print 'Time:\t\t\t {}'.format(time.clock() - start_total)
    print 'Output file:\t {}'.format(OUTPUT_FILE)

    print_result(OUTPUT_FILE)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("f", type=str, help="path to input file with data")
    parser.add_argument("X", type=int, help="width of buffer zone rectangle")
    parser.add_argument("Y", type=int, help="height of buffer zone rectangle")
    parser.add_argument("az", type=int, help="azimuth threshold")
    parser.add_argument("o", type=str, help="output file")
    args = parser.parse_args()

    VERTICES = read_data(args.f)
    VERTICES_LEN = len(VERTICES)
    BORDER_X = args.X
    BORDER_Y = args.Y
    BORDER_AZ = args.az
    OUTPUT_FILE = args.o

    clustering()
