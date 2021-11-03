/*******************************************************************************
//  3TSPinsertion.cpp               Author: Ian Nobile
//  Section: 50                     Due Date: 24 September 2021
//
//  This memory-leak- and bug-free program solves provided instances of a 
//  Travelling Salesperson Problem by using a Closest Edge Insertion Heuristic, 
//  a variant of the greedy heuristic. It builds its final tour by starting with 
//  a single city and expanding to three. It then cycles through each edge of the
//  graph, checking distances to each of the remaining unvisited points before 
//  inserting those corresponding to the shortest distances one by one.
//
*******************************************************************************/

#include <iostream> // print to console
#include <chrono>   // time the speed of the program
#include <fstream>  // read files
#include <string>   // header buffer for seeking inside a file
#include <cfloat>   // allows the use of FLT_MAX
#include <vector>   // easy arraying
#include <cmath>    // distance formula
#include <C:\Users\admin\Desktop\coastal_islander\me\uofl\ai\project3\reso\matplotlib-cpp-master\matplotlibcpp.h>

using namespace std;
using namespace std::chrono;
namespace plt = matplotlibcpp;


// class declarations:
class Node {
public:
    int num;
    float x;
    float y;
    float dist;     // ?
    Node *child;
    bool visited;   // ?
};

class Graph {
public:
    vector<Node> nodes;
};


// function prototypes:
Graph buildGraph(char *);
float distCheck(Node, Node);
float distCheck(Node, Node, Node);
void tourAdopter(vector<Node> &);



//------------------------------------------------------------------------------
//  Main Function
//------------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    // check if path was passed as arg:
    if (argc == 1) {
        cout << "Please pass the path to the .TSP file as a command line argument" << endl;
        return 1;
    }

    // begin with a friendly greeting:
    cout << "Hello and welcome to the (Closest Edge Insertion Heuristic) Travelling Salesperson Problem Solver" << endl;

    //Fill graph
    Graph graph = buildGraph(argv[1]);
    vector<Node> tour;

    // Graph the points heere!
    /*
    vector<double> xPlot;
    vector<double> yPlot;
    plt::xkcd();
    for (Node node : graph.nodes) {
        xPlot.push_back(node.x);
        yPlot.push_back(node.y);
        plt::plot(xPlot, yPlot,"-o");
        xPlot.pop_back();
        yPlot.pop_back();
    }
    plt::show();
    */
    auto start = high_resolution_clock::now();  // start timer

    //Insert 0 into tour
    vector<Node>::iterator graphIt = graph.nodes.begin();

    tour.push_back(*graphIt);
    graph.nodes.erase(graphIt);
    float minDist = FLT_MAX;
    float currDist = FLT_MAX;
    vector<Node>::iterator tourIt = tour.begin();
    Node closest = Node();

    //Check all 2-D dists between all points and 0
    for (graphIt = graph.nodes.begin();graphIt != graph.nodes.end();graphIt++) {
        currDist = distCheck(*tourIt, *graphIt);
        if (currDist < minDist) {
            minDist = currDist;
            closest = *graphIt;
        }
    }
    graphIt = graph.nodes.begin();
    while (graphIt->num != closest.num) { graphIt++; }

    //Insert shortest into tour
    tour.push_back(*graphIt);
    tourAdopter(tour);
    graph.nodes.erase(graphIt);
    minDist = FLT_MAX;
    currDist = FLT_MAX;
    tourIt = tour.begin();

    //Check all 3-D dists between all points
    for (graphIt = graph.nodes.begin();graphIt != graph.nodes.end();graphIt++) {
        currDist = distCheck(*tourIt, *tourIt->child, *graphIt);
        if (currDist < minDist) {
            minDist = currDist;
            closest = *graphIt;
        }
    }
    graphIt = graph.nodes.begin();
    while (graphIt->num != closest.num) { graphIt++; }

    //Insert shortest into tour
    tour.push_back(*graphIt);
    tourAdopter(tour);
    graph.nodes.erase(graphIt);
    Node parent = Node();

    // cycle through all tour edges checking 3-D dists between all points
    while (graph.nodes.size() > 0) {
        minDist = FLT_MAX;
        currDist = FLT_MAX;
        for (tourIt = tour.begin();tourIt != tour.end();tourIt++) {
            for (graphIt = graph.nodes.begin();graphIt != graph.nodes.end();graphIt++) {
                currDist = distCheck(*tourIt, *tourIt->child, *graphIt);
                if (currDist < minDist) {
                    minDist = currDist;
                    closest = *graphIt;
                    parent = *tourIt;
                }
            }
        }
        graphIt = graph.nodes.begin();
        while (graphIt->num != closest.num) { graphIt++; }
        tourIt = tour.begin();
        while (tourIt->num != parent.child->num) { tourIt++; }

        //Insert shortest between parent and child in tour;
        if (tourIt->num == tour.begin()->num) {
            tour.push_back(*graphIt);
        } else {
            tour.insert(tourIt, *graphIt);
        }
        tourAdopter(tour);
        graph.nodes.erase(graphIt);
        // repeat until all edges added to tour and graph is empty
    }
    
    // Graph the points and calculate minDist heere!
    cout << "The shortest route visiting all the points is:" << endl << endl;
    minDist = 0;
    for (int i = 0;i < tour.size() - 1;i++) {
        cout << tour[i].num;
        //xPlot.push_back(tour[i].x);
        //yPlot.push_back(tour[i].y);
        //plt::plot(xPlot, yPlot);
        //plt::pause(1);
        minDist += distCheck(tour[i], tour[i + 1]);
        if (i == tour.size() - 1) { break; }
        cout << " -> ";
    }
    cout << tour[tour.size() - 1].num;
    //xPlot.push_back(tour[tour.size() - 1].x);
    //yPlot.push_back(tour[tour.size() - 1].y);
    //plt::plot(xPlot, yPlot);
    //plt::pause(1);
    minDist += distCheck(tour[tour.size() - 1], tour[0]);
    //xPlot.push_back(tour[0].x);
    //yPlot.push_back(tour[0].y);
    //plt::plot(xPlot, yPlot);
    //plt::pause(1);
    //plt::show();
    cout << endl << endl << "and back again at a total distance of " << minDist << endl << endl;

    auto stop = high_resolution_clock::now();   // stop timer
    auto duration = duration_cast<microseconds>(stop - start);  // calculate elapsed time
    cout << "Program execution took " << duration.count() / 1000000.0 << "s" << endl << endl;

    system("pause");
    return 0;
}


//  function definitions:
//------------------------------------------------------------------------------
//  Read .TSP file, creates nodes from coordinates and combines all in an
//  undirected graph object
//------------------------------------------------------------------------------
Graph buildGraph(char *argv) {
    // open .TSP file in read-only mode:
    ifstream tspfile;
    tspfile.open(argv, ios::in);
    // ensure file exists:
    if (!tspfile.is_open()) {
        Graph graph;
        return graph;
    }
    Graph graph;
    int dimension = 0;
    string heading = "";
    // advance buffer to the dimension section:
    while (heading.compare("DIMENSION:") != 0) {
        tspfile >> heading;
    }
    tspfile >> dimension;
    // advance buffer to the coordinates section:
    while (heading.compare("NODE_COORD_SECTION") != 0) {
        tspfile >> heading;
    }
    // create nodes and push to graph vector
    Node newNode = Node();
    for (int i = 0;i < dimension;i++) {
        tspfile >> newNode.num;
        tspfile >> newNode.x;
        tspfile >> newNode.y;
        newNode.dist = FLT_MAX;
        // newNode.child;
        newNode.visited = false;
        graph.nodes.push_back(newNode);
    }
    // The graph is now created, and we are finished with the .TSP file
    tspfile.close();
    return graph;
}

//------------------------------------------------------------------------------
//  Returns the distance between two graph nodes using the pythagoran formula: 
//  dist = sqrt((x2 - x1)^2 + (y2 - y1)^2)
//------------------------------------------------------------------------------
float distCheck(Node first, Node second) {
    return sqrt((pow((second.x - first.x), 2.0)) + (pow((second.y - first.y), 2.0)));
}

//------------------------------------------------------------------------------
//  Returns the distance between three graph nodes by creating vectors and
//  choosing an appropriate distance calculation method
//------------------------------------------------------------------------------
float distCheck(Node parent, Node child, Node point) {
    
    //adapted from pDistance() by Joshua Perina:
    float A = point.x - parent.x;
    float B = point.y - parent.y;
    float C = child.x - parent.x;
    float D = child.y - parent.y;

    float dot = A * C + B * D;
    float len_sq = C * C + D * D;
    int param = dot / len_sq;

    float xx, yy;

    if (param < 0) {
        xx = parent.x;
        yy = parent.y;
    }
    else if (param > 1) {
        xx = child.x;
        yy = child.y;
    }
    else {
        xx = parent.x + param * C;
        yy = parent.y + param * D;
    }
    float dx = point.x - xx;
    float dy = point.y - yy;
    return sqrt(dx * dx + dy * dy);
}

//------------------------------------------------------------------------------
//  Returns the distance between three graph nodes using the pythagoran formula: 
//  dist = sqrt((x2 - x1)^2 + (y2 - y1)^2)
//------------------------------------------------------------------------------
void tourAdopter(vector<Node> &tour) {
    for (int i = 0;i < tour.size() - 1;i++) {
        tour[i].child = &tour[i + 1];
    }
    tour[tour.size() - 1].child = &tour[0];
    return;
}

