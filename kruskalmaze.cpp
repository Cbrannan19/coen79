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