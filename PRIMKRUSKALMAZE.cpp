#include <iostream>
#include <stack>
#include "maze.hpp"


// Construtor for Maze Class
// Input: Int, Int representing row and columns in maze
// Outout: Maze Object of size r x c

Maze::Maze(int r, int c) {
    if (r % 2 == 0)
        r++;
    if (c % 2 == 0)
        c++;
    row = r;
    col = c;
    int paddedRows = (r * 2) + 1;
    int paddedCols = (c * 2) + 1;

    std::vector<std::tuple<int, char>> rows;
    int idx = 0;
    for (int i = 0; i < paddedRows; ++i) {
        for (int j = 0; j < paddedCols; ++j) {
            if (i % 2 == 1 && j % 2 == 1) {
                auto t = std::make_tuple(idx, '0');
                idx++;
                rows.push_back(t);
            } else {
                auto t = std::make_tuple(-1, '1');
                rows.push_back(t);
            }
        }
        grid.push_back(rows);
        rows.erase(rows.begin(), rows.end());
    }
    numRows = paddedRows;
    numCols = paddedCols;
    numVertices = row * col;
    makeEdges();
    kGrid = grid;
    pGrid = grid;
    for(int i = 0; i < grid.size(); ++i){
        std::vector<bool> temp (numCols, false);
        kPathGrid.push_back(temp);
        pPathGrid.push_back(temp);
    }
}


// This is my creative way to generate weights for this. Instead of using just rand()
// or srand() I am using a PRNG to generate a random number between 1 - 100 based off the srand() and rand() seed
// passed to it. To ensure more variation, it is rerun every time a weight is needed so the pool of available numbers
// doesn't decrease.
// Input: N/A
// Output: N/A

void Maze::makeEdges() {
    std::srand(time(NULL));
    int vertexFrom = 0;

    for (int i = 0; i < getNumRows(); ++i) {
        for (int j = 0; j < getNumCols(); ++j) {
            if (i % 2 == 1 && j % 2 == 1 && i != 0 && i != numRows - 1 && j != 0 && j != numCols - 1) {
                if (j < getNumCols() - 2) {
                    int vertexTo = vertexFrom + 1;
                    int weight = genRandWeight();
                    auto fromh = std::make_tuple(i, j);
                    auto toh = std::make_tuple(i, j + 2);
                    auto edgeh = std::make_tuple(fromh, toh, weight, vertexFrom, vertexTo, "horiz");
                    edges.push_back(edgeh);
                }
                if (i < getNumRows() - 2) {
                    int vertexTo = vertexFrom + col;
                    int weight = genRandWeight();
                    auto fromv = std::make_tuple(i, j);
                    auto tov = std::make_tuple(i + 2, j);
                    auto edgev = std::make_tuple(fromv, tov, weight, vertexFrom, vertexTo, "vert");
                    edges.push_back(edgev);
                }
                ++vertexFrom;
            }
        }
    }
    source();
    dest();
    while (std::get<3>(sourceEdge) == std::get<3>(destEdge) || std::get<4>(sourceEdge) == std::get<4>(destEdge)) {
        source();
    }
}


// Generates a random weight for the maze edges
// Input: N/A
// Output: A random Int from A Pseudo Random Number Generator

int Maze::genRandWeight() {
    int seed = rand() % 16;
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> distribution(1, 100);
    int weight = distribution(generator);
    return weight;
}


// Randomly picks source of maze using a PRNG. Check for modulus and that determines
// what edge to start on
// Input: N/A
// Output: tuple with the row, column location of the start point

std::tuple<int, int> Maze::genRandSourcePos() {
    int seed = rand() % 100;
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> distribution(1, 1000);
    int sideOrTop = distribution(generator) % 2;
    if (sideOrTop == 0) {
        int row = distribution(generator);
        row = row % (numRows - 2);
        if (row % 2 == 0)
            row++;
        auto pos = std::make_tuple(row, 0);
        std::get<1>(grid.at(row).at(0)) = '0';
        return pos;
    } else {
        int col = distribution(generator);
        col = col % (numCols - 2);
        if (col % 2 == 0)
            col++;
        auto pos = std::make_tuple(0, col);
        std::get<1>(grid.at(0).at(col)) = '0';
        return pos;
    }
}


// Randomly picks source of maze using a PRNG. Check for modulus and that determines
// what edge to end on
// Input: N/A
// Output: tuple with the row, column location of the end point

std::tuple<int, int> Maze::genRandDestPos() {
    int seed = rand() % 100;
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> distribution(1, 1000);
    int sideOrTop = distribution(generator) % 2;
    if (sideOrTop == 0) {
        int row = distribution(generator);
        row = row % (numRows - 2);
        auto t = std::get<0>(sourceEdge);
        while (row == std::get<0>(t)) {
            row = distribution(generator);
            row = row % (numRows - 2);
        }
        if (row % 2 == 0)
            row++;
        auto pos = std::make_tuple(row, numCols - 1);
        std::get<1>(grid.at(row).at(numCols - 1)) = '0';
        return pos;
    } else {
        int col = distribution(generator);
        col = col % (numCols - 2);
        auto t = std::get<1>(sourceEdge);
        while (col == std::get<1>(t)) {
            col = distribution(generator);
            col = col % (numCols - 2);
        }
        if (col % 2 == 0)
            col++;
        auto pos = std::make_tuple(numRows - 1, col);
        std::get<1>(grid.at(numRows - 1).at(col)) = '0';
        return pos;
    }
}


// Creates start position and sets weight to -10000 to keep it from being mistaken for another weighted edge
// Input: N/A
// Output: N/A

void Maze::source() {
    int weight = genRandWeight();
    auto from = genRandSourcePos();
    if (std::get<0>(from) == 0) {
        auto to = std::make_tuple(std::get<0>(from) + 1, std::get<1>(from));
        int x = (std::get<0>(grid.at(std::get<0>(to)).at(std::get<1>(to))));
        sourceEdge = std::make_tuple(from, to, weight, -10000, x, "source");
    } else {
        auto to = std::make_tuple(std::get<0>(from), std::get<1>(from) + 1);
        int x = (std::get<0>(grid.at(std::get<0>(to)).at(std::get<1>(to))));
        sourceEdge = std::make_tuple(from, to, weight, -10000, x, "source");
    }
}


// Creates end position and sets weight to -10000 to keep it from being mistaken for another weighted edge
// Input: N/A
// Output: N/A

void Maze::dest() {
    int weight = genRandWeight();
    auto to = genRandDestPos();
    if (std::get<0>(to) == numRows - 1) {
        auto from = std::make_tuple(std::get<0>(to) - 1, std::get<1>(to));
        int x = (std::get<0>(grid.at(std::get<0>(from)).at(std::get<1>(from))));
        destEdge = std::make_tuple(from, to, weight, x, 10000, "dest");
    } else {
        auto from = std::make_tuple(std::get<0>(to), std::get<1>(to) - 1);
        int x = (std::get<0>(grid.at(std::get<0>(from)).at(std::get<1>(from))));
        destEdge = std::make_tuple(from, to, weight, x, 10000, "dest");
    }
}


// Calls the Implementation of Prim's Algorithm for a minimum spanning tree to generate the edges in the maze
// Input: N/A
// Output: N/A

void Maze::prims() {
    std::vector<bool> visited(numVertices, false);
    std::vector<int> v;
    std::vector<std::vector<int>> adjList;
    pEdgesFinal.push_back(sourceEdge);
    primsMST(adjList, visited);
    pEdgesFinal.push_back(destEdge);
    std::cout << "Prims's Maze:\n\n";
    print(pGrid);
    auto t = std::get<0>(sourceEdge);
    pShortestPath(std::get<0>(t), std::get<1>(t));
    std::cout << "Prim's Maze Solved:\n\n";
    print(pGrid);
}


// Implementation of Prim's Algorithm for a minimum spanning tree to generate the edges in the maze
// Input: N/A
// Output: N/A


void Maze::primsMST(std::vector<std::vector<int>> adjList, std::vector<bool> visited) {
    std::vector<std::vector<std::tuple<std::tuple<int, int>, std::tuple<int, int>, int, int, int, std::string>>> heap;
    std::vector<std::vector<std::tuple<std::tuple<int, int>, std::tuple<int, int>, int, int, int, std::string>>> temp;

    // Populate a temp vector of edges and repopulate the queue
    while (!primsEdges.empty()) {
        std::vector<std::tuple<std::tuple<int, int>, std::tuple<int, int>, int, int, int, std::string>> t;
        t.push_back(primsEdges.front());
        temp.push_back(t);
        primsEdges.pop();
    }
    for (int i = 0; i < temp.size(); ++i) {
        for (int j = 0; j < temp.at(i).size(); ++j) {
            primsEdges.push(temp.at(i).at(j));
        }
    }
    int pos = std::get<4>(sourceEdge);
    visited.at(pos) = true;
    while (!areAllVisited(visited)) {
        for (int i = 0; i < temp.size(); ++i) {
            for (int j = 0; j < temp.at(i).size(); ++j) {
                if (std::get<4>(temp.at(i).at(j)) == pos || std::get<3>(temp.at(i).at(j)) == pos) {
                    heap.push_back(temp.at(i));
                    std::get<3>(temp.at(i).at(j)) = 999;
                    std::get<4>(temp.at(i).at(j)) = 999;
                }
            }
        }
        std::tuple<std::tuple<int, int>, std::tuple<int, int>, int, int, int, std::string> minWeightedEdge;//FIXME: run through to fine tune
        std::get<2>(minWeightedEdge) = 101;
        for (int i = 0; i < heap.size(); ++i) {
            for (int j = 0; j < heap.at(i).size(); ++j) {
                if (std::get<2>(minWeightedEdge) > std::get<2>(heap.at(i).at(j))) {
                    minWeightedEdge = heap.at(i).at(j);
                }
            }
        }
        if (!hasCycle(std::get<3>(minWeightedEdge), std::get<4>(minWeightedEdge), adjList, visited)) {
            pEdgesFinal.push_back(minWeightedEdge);
            if (std::get<5>(minWeightedEdge) == "horiz") {
                auto t = std::get<0>(minWeightedEdge);
                std::get<1>(pGrid.at(std::get<0>(t)).at(std::get<1>(t) + 1)) = '0';
            }
            if (std::get<5>(minWeightedEdge) == "vert") {
                auto t = std::get<0>(minWeightedEdge);
                std::get<1>(pGrid.at(std::get<0>(t) + 1).at(std::get<1>(t))) = '0';
            }
            visited.at(std::get<3>(minWeightedEdge)) = true;
            visited.at(std::get<4>(minWeightedEdge)) = true;
            std::vector<int> v;
            v.push_back(pos);
            if (adjList.empty()) {
                adjList.push_back(v);
                v.clear();
                if(std::get<4>(minWeightedEdge) == pos){
                    adjList.at(0).push_back(std::get<3>(minWeightedEdge));
                }
                else {
                    adjList.at(0).push_back(std::get<4>(minWeightedEdge));
                }
            }
            else {
                int fromIdx = -9, toIdx = -8;
                for (int i = 0; i < adjList.size(); ++i) {
                    for (int j = 0; j < adjList.at(i).size(); ++j) {
                        if (std::get<4>(minWeightedEdge) == adjList.at(i).at(j)) {
                            fromIdx = i;
                        }
                        if (std::get<3>(minWeightedEdge) == adjList.at(i).at(j)) {
                            toIdx = i;
                        }
                    }
                }
                if (fromIdx != -9 && toIdx == -8) {
                    bool existsTo = false, existsFrom = false;
                    for (int i = 0; i < adjList.at(fromIdx).size(); ++i) {
                        if (adjList.at(fromIdx).at(i) == std::get<3>(minWeightedEdge))
                            existsFrom = true;
                        if (adjList.at(fromIdx).at(i) == std::get<4>(minWeightedEdge))
                            existsTo = true;
                    }
                    if(!existsTo) {
                        adjList.at(fromIdx).push_back(std::get<4>(minWeightedEdge));
                    }
                    if(!existsFrom) {
                        adjList.at(fromIdx).push_back(std::get<3>(minWeightedEdge));
                    }
                } else if (fromIdx == -9 && toIdx != -8) {
                    bool existsTo = false, existsFrom = false;
                    for (int i = 0; i < adjList.at(toIdx).size(); ++i) {
                        if (adjList.at(toIdx).at(i) == std::get<3>(minWeightedEdge))
                            existsFrom = true;
                        if (adjList.at(toIdx).at(i) == std::get<4>(minWeightedEdge))
                            existsTo = true;
                    }
                    if(!existsTo) {
                        adjList.at(toIdx).push_back(std::get<4>(minWeightedEdge));
                    }
                    if(!existsFrom) {
                        adjList.at(toIdx).push_back(std::get<3>(minWeightedEdge));
                    }
                }
            }
            pos = adjList.at(0).back();
            auto iter = temp.begin();
            for(int i = 0; i < temp.size(); ++i){
                for(int j = 0; j < temp.at(i).size(); ++j){
                    if(std::get<3>(minWeightedEdge) == std::get<3>(temp.at(i).at(j)) &&
                            std::get<4>(minWeightedEdge) == std::get<4>(temp.at(i).at(j))){
                        temp.erase(iter);
                        break;
                    }
                }
                ++iter;
            }
            iter = heap.begin();
            for(int i = 0; i < heap.size(); ++i) {
                for (int j = 0; j < heap.at(i).size(); ++j) {
                    if(std::get<3>(minWeightedEdge) == std::get<3>(heap.at(i).at(j)) &&
                    std::get<4>(minWeightedEdge) == std::get<4>(heap.at(i).at(j))){
                        heap.erase(iter);
                        break;
                    }
                }
                ++iter;
            }
        }
        else{
            auto iter = temp.begin();
            for(int i = 0; i < temp.size(); ++i){
                for(int j = 0; j < temp.at(i).size(); ++j){
                    if(std::get<3>(minWeightedEdge) == std::get<3>(temp.at(i).at(j)) &&
                       std::get<4>(minWeightedEdge) == std::get<4>(temp.at(i).at(j))){
                        temp.erase(iter);
                        break;
                    }
                }
                ++iter;
            }
            iter = heap.begin();
            for(int i = 0; i < heap.size(); ++i) {
                for (int j = 0; j < heap.at(i).size(); ++j) {
                    if(std::get<3>(minWeightedEdge) == std::get<3>(heap.at(i).at(j)) &&
                       std::get<4>(minWeightedEdge) == std::get<4>(heap.at(i).at(j))){
                        heap.erase(iter);
                        break;
                    }
                }
                ++iter;
            }
        }
    }
}


// Keeps tabs if all edges have been visited
// Input: vector of booleans to keep track of visited edges

bool Maze::areAllVisited(std::vector<bool> visited) {
    for (int i = 0; i < visited.size(); ++i) {
        if (!visited.at(i))
            return false;
    }
    return true;
}


// Calls the Implementation of Kruskal's Algorithm for a minimum spanning tree to generate the edges in the maze
// Input: N/A
// Output: N/A

void Maze::kruskals() {
    std::vector<bool> visited(numVertices, false);
    std::vector<int> v;
    std::vector<std::vector<int>> adjList;
    kruskalSort();
    primsEdges = kruskalsEdges;
    kEdgesFinal.push_back(sourceEdge);
    kruskalsMST(adjList, visited);
    kEdgesFinal.push_back(destEdge);
    std::cout << "Kruskal's Maze:\n\n";
    print(kGrid);
    auto t = std::get<0>(sourceEdge);
    kShortestPath(std::get<0>(t), std::get<1>(t));
    std::cout << "Kruskal's Maze Solved:\n\n";
    print(kGrid);
}


// Sorting algorithm to sort the kruskal vector of edges by weight
// Input: N/A
// Output: N/A

void Maze::kruskalSort() {
    std::vector<std::tuple<std::tuple<int, int>, std::tuple<int, int>, int, int, int, std::string>> temp = edges;
    std::stack<std::tuple<std::tuple<int, int>, std::tuple<int, int>, int, int, int, std::string>> stack;
    while (stack.size() < edges.size()) {
        int idx = -1;
        int min = 0;
        for (int i = 0; i < temp.size(); ++i) {
            if (std::get<2>(temp.at(i)) > min) {
                min = std::get<2>(temp.at(i));
                idx = i;
            }
        }
        stack.push(temp.at(idx));
        temp.erase(temp.begin() + idx);
    }
    while (stack.size() > 0) {
        auto t = stack.top();
        stack.pop();
        kruskalsEdges.push(t);
    }
}

// Implementation of Kruskal's Algorithm for a minimum spanning tree to generate the edges in the maze
// Input: N/A
// Output: N/A



void Maze::kruskalsMST(std::vector<std::vector<int>> adjList, std::vector<bool> visited) {
    while (!kruskalsEdges.empty()) {
        auto e = kruskalsEdges.front();
        if (!hasCycle(std::get<3>(e), std::get<4>(e), adjList, visited)) {
            kEdgesFinal.push_back(kruskalsEdges.front());
            if (std::get<5>(e) == "horiz") {
                auto t = std::get<0>(e);
                std::get<1>(kGrid.at(std::get<0>(t)).at(std::get<1>(t) + 1)) = '0';
            }
            if (std::get<5>(e) == "vert") {
                auto t = std::get<0>(e);
                std::get<1>(kGrid.at(std::get<0>(t) + 1).at(std::get<1>(t))) = '0';
            }
            std::vector<int> temp;
            temp.push_back(std::get<3>(e));
            if (adjList.empty()) {
                adjList.push_back(temp);
                adjList.at(0).push_back(std::get<4>(e));
            } else {
                int fromIdx = -9, toIdx = -8;
                for (int i = 0; i < adjList.size(); ++i) {
                    for (int j = 0; j < adjList.at(i).size(); ++j) {
                        if (adjList.at(i).at(j) == std::get<3>(e)) {
                            fromIdx = i;
                        }
                        if (adjList.at(i).at(j) == std::get<4>(e)) {
                            toIdx = i;
                        }
                    }
                }
                if (fromIdx == toIdx) {
                    break;
                } else if (fromIdx != -9 && toIdx == -8) {
                    adjList.at(fromIdx).push_back(std::get<4>(e));
                } else if (fromIdx == -9 && toIdx != -8) {
                    adjList.at(toIdx).push_back(std::get<3>(e));
                } else if (fromIdx != -9 && toIdx != -8) {
                    if (adjList.at(fromIdx).size() < adjList.at(toIdx).size()) {
                        while (!adjList.at(fromIdx).empty()) {
                            int x = adjList.at(fromIdx).back();
                            adjList.at(fromIdx).pop_back();
                            adjList.at(toIdx).push_back(x);
                        }
                        adjList.erase(adjList.begin() + fromIdx);
                    } else {
                        while (!adjList.at(toIdx).empty()) {
                            int x = adjList.at(toIdx).back();
                            adjList.at(toIdx).pop_back();
                            adjList.at(fromIdx).push_back(x);
                        }
                        adjList.erase(adjList.begin() + toIdx);
                    }
                }
                if (fromIdx == -9 && toIdx == -8) {
                    temp.push_back(std::get<4>(e));
                    adjList.push_back(temp);
                }
            }
            visited.at(std::get<3>(e)) = visited.at(std::get<4>(e)) = true;
        }
        kruskalsEdges.pop();
    }
}


// Checks for cycles in the maze
// Input: Int, Int, vector of Ints of adjacent edges and a vector of Booleans for visited vertices
// Output: Boolean which signifies whether a cycle exists(true) or not(false)

bool Maze::hasCycle(int from, int to, std::vector<std::vector<int>> adjList, std::vector<bool> visited) {
    if (visited.at(to)) {
        int fromIdx = -9, toIdx = -8;
        for (int i = 0; i < adjList.size(); ++i) {
            for (int j = 0; j < adjList.at(i).size(); ++j) {
                if (adjList.at(i).at(j) == from)
                    fromIdx = i;
                if (adjList.at(i).at(j) == to)
                    toIdx = i;
                if (fromIdx == toIdx)
                    return true;
            }
        }
    }
    return false;
}


// Finds the shortest path of the Prim's Algorithm maze
// Input: Int, Int: row and column of the current vertex
// Output: Boolean which signifies if the current maze was solved or not

bool Maze::pShortestPath(int rowx, int colx) {
    auto t = std::get<1>(destEdge);
    if(rowx < 0 || colx < 0) {
        return false;
    }
    if(rowx >= numRows || colx >= numCols) {
        return false;
    }
    if(rowx == std::get<0>(t) && colx == std::get<1>(t)) {
        std::get<1>(pGrid.at(rowx).at(colx)) = ' ';
        return true;
    }
    if(std::get<1>(pGrid.at(rowx).at(colx)) == '1'){
        return false;
    }
    if(pPathGrid.at(rowx).at(colx)) {
        return false;
    }

    pPathGrid.at(rowx).at(colx) = true;

    if(pShortestPath(rowx, colx + 1)){
        std::get<1>(pGrid.at(rowx).at(colx)) = ' ';
        return true;
    }
    if(pShortestPath(rowx, colx - 1)){
        std::get<1>(pGrid.at(rowx).at(colx)) = ' ';
        return true;
    }
    if(pShortestPath(rowx - 1, colx)){
        std::get<1>(pGrid.at(rowx).at(colx)) = ' ';
        return true;
    }
    if(pShortestPath(rowx + 1, colx)){
        std::get<1>(pGrid.at(rowx).at(colx)) = ' ';
        return true;
    }

    pPathGrid.at(rowx).at(colx) = false;
    return false;
}


// Finds the shortest path of the Kruskal's Algorithm maze
// Input: Int, Int: row and column of the current vertex
// Output: Boolean which signifies if the current maze was solved or not

bool Maze::kShortestPath(int rowx, int colx) {
    auto t = std::get<1>(destEdge);
    if(rowx < 0 || colx < 0) {
        return false;
    }
    if(rowx >= numRows || colx >= numCols) {
        return false;
    }
    if(rowx == (std::get<0>(t)) && colx == std::get<1>(t)) {
        std::get<1>(kGrid.at(rowx).at(colx)) = ' ';
        return true;
    }
    if(std::get<1>(kGrid.at(rowx).at(colx)) == '1'){
        return false;
    }
    if(kPathGrid.at(rowx).at(colx)) {
        return false;
    }
    kPathGrid.at(rowx).at(colx) = true;

    if(kShortestPath(rowx, colx + 1)){
        std::get<1>(kGrid.at(rowx).at(colx)) = ' ';
        return true;
    }
    if(kShortestPath(rowx, colx - 1)){
        std::get<1>(kGrid.at(rowx).at(colx)) = ' ';
        return true;
    }
    if(kShortestPath(rowx - 1, colx)){
        std::get<1>(kGrid.at(rowx).at(colx)) = ' ';
        return true;
    }
    if(kShortestPath(rowx + 1, colx)){
        std::get<1>(kGrid.at(rowx).at(colx)) = ' ';
        return true;
    }

    kPathGrid.at(rowx).at(colx) = false;
    return false;
}


// Prints the maze to the console
// Input: a vector of vectors of tuples of Int, Char that is the maze to be printed
// Output: N/A

void Maze::print(std::vector<std::vector<std::tuple<int, char>>> grid) {
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            std::cout << std::get<1>(grid.at(i).at(j)) << " ";
            if (j == numCols - 1)
                std::cout << std::endl;
        }
    }
    std::cout << std::endl << std::endl;
}