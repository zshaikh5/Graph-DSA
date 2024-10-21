// Jasleen Kaur Saini and Zaina Shaikh 
// CSS343 Homework 3: Graph 

#include "graph.h"
#include <algorithm>
#include <climits>
#include <fstream>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>
#include <functional>
#include <iostream>


//using namespace std;

////////////////////////////////////////////////////// constructor + destructor ////////////////////////////////////////////////////////////

// constructor, empty graph
// directionalEdges defaults to true
Graph::Graph(bool directionalEdges) 
{ 
  directional = !directionalEdges; 
}

// destructor
Graph::~Graph() 
{
    for (int i = 0; i < vertexList.size(); i++) 
    {
        delete vertexList[i];
    }
    vertexList.clear();
}

/////////////////////////////////////////////////////////////// functions ///////////////////////////////////////////////////////////////

// read a text file and create the graph
bool Graph::readFile(const string &filename) 
{
  ifstream myfile(filename);
  if (!myfile.is_open()) 
  {
    cerr << "Failed to open " << filename << endl;
    return false;
  }
  int edges = 0;
  int weight = 0;
  string fromVertex;
  string toVertex;
  myfile >> edges;
  for (int i = 0; i < edges; ++i) 
  {
    myfile >> fromVertex >> toVertex >> weight;
    connect(fromVertex, toVertex, weight);
  }
  myfile.close();
  return true;
}

// @return total number of vertices
int Graph::verticesSize() const 
{ 
  return vertexList.size(); 
}

// @return total number of edges
int Graph::edgesSize() const 
{
  int temp = 0;
  for (int i=0; i<vertexList.size(); i++) 
  {
    temp += vertexList[i]->edgeList.size();
  }
  if (directional) 
  {
    temp /= 2;
  }
  return temp;
}

// @return integer index value of given vertex 
int Graph::index(const string &from, const string &label) const 
{
    if (from.empty()) 
    {
        for (int i = 0; i < vertexList.size(); i++) 
        {
            if (vertexList[i]->value == label) 
            {
              return i; 
            }
        }
        return -1;
    }

    int inVertexIndex = index("", from);
    if (inVertexIndex == -1) 
    {
      return -1;
    }

    for (int i = 0; i < vertexList[inVertexIndex]->edgeList.size(); i++) 
    {
      if (vertexList[inVertexIndex]->edgeList[i].first->value == label) 
      {
        return i; 
      }
    }
    return -1; 
}

// @return number of edges from given vertex, -1 if vertex not found
int Graph::vertexDegree(const string &label) const 
{
    int vertexIndex = index("", label);
    if (vertexIndex != -1) 
    {
        int degree = vertexList[vertexIndex]->edgeList.size();
        return degree;
    }
    return vertexIndex;
}

// @return true if vertex added, false if it already is in the graph
bool Graph::add(const string &label) 
{
    if (!contains(label)) 
    {
        Vertex *newVertex = new Vertex(label);
        bool inserted = false;

        for (int i = 0; i < vertexList.size(); i++) 
        {
            if (vertexList[i]->value > label) 
            {
                vertexList.insert(vertexList.begin() + i, newVertex);
                inserted = true;
                break; 
            }
        }

        if (!inserted) 
        {
            vertexList.push_back(newVertex);
        }

        return true; 
    }
    return false;
}

/** return true if vertex already in graph */
bool Graph::contains(const string &label) const 
{
    int indexResult = index("", label);
    if (indexResult != -1) 
    {
      return true; 
    } 
    else 
    {
        return false; 
    }
}

// @return string representing edges and weights, "" if vertex not found
// A-3->B, A-5->C should return B(3),C(5)
string Graph::getEdgesAsString(const string &label) const 
{
    int indexNum = index("", label);
    if (indexNum != -1 && !vertexList[indexNum]->edgeList.empty()) 
    {
        string edge;
        vector<pair<Vertex*, int>>& edges = vertexList[indexNum]->edgeList;
        for (int i = 0; i < edges.size(); ++i) 
        {
            edge.append(edges[i].first->value + "(" + to_string(edges[i].second) + "),");
        }
        edge.erase(edge.end() - 1);
        return edge;
    }
    return "";
}

// @return true if successfully connected
bool Graph::connect(const string &from, const string &to, int weight) 
{
    if (!contains(from)) add(from);
    if (!contains(to) && from != to) add(to);

    int fromIndex = index("", from);
    int toIndex = index("", to);

    if (index(from, to) != -1 || from == to) 
    {
      return false;
    }

    vector<pair<Vertex *, int>>& fromEdges = vertexList[fromIndex]->edgeList;
    vector<pair<Vertex *, int>>::iterator it = 
      lower_bound(fromEdges.begin(), fromEdges.end(), to,
      [](const pair<Vertex *, int>& lhs, const string& rhs) 
      { return lhs.first->value < rhs;});

    fromEdges.insert(it, make_pair(vertexList[toIndex], weight));

    if (directional) 
    {
        vector<pair<Vertex *, int>>& toEdges = vertexList[toIndex]->edgeList;
        vector<pair<Vertex *, int>>::iterator it = 
          lower_bound(toEdges.begin(), toEdges.end(), from,
          [](const pair<Vertex *, int>& lhs, const string& rhs) {
          return lhs.first->value < rhs;});
        toEdges.insert(it, make_pair(vertexList[fromIndex], weight));
    }
    return true;
}

bool Graph::disconnect(const string &from, const string &to) 
{
    if (!contains(from) || !contains(to) || from == to) 
    {
        return false;
    }

    int fromIndex = index("", from);
    int toIndex = index("", to);

    vector<pair<Vertex*, int>>& fromEdges = vertexList[fromIndex]->edgeList;
    for (int i = 0; i < fromEdges.size(); ++i) 
    {
        if (fromEdges[i].first->value == to) 
        {
            fromEdges.erase(fromEdges.begin() + i);
            if (directional) 
            {
                vector<pair<Vertex*, int>>& toEdges = vertexList[toIndex]->edgeList;
                for (int j = 0; j < toEdges.size(); ++j) 
                {
                    if (toEdges[j].first->value == from) 
                    {
                        toEdges.erase(toEdges.begin() + j);
                        break;
                    }
                }
            }
            return true;
        }
    }
    return false;
}

//////////////////////////////////////////////////////////// algorithms ///////////////////////////////////////////////////////////////

// depth-first traversal starting from given startLabel
void Graph::dfs(const string &startLabel, void (*visit)(const string &label)) 
{
    int startIndex = index("", startLabel);
    if (startIndex == -1) return; 
    vector<string> visited; 
    traverseVerticesHelper(vertexList[startIndex], visited); 
    for (int i = 0; i < visited.size(); i++) 
    {
        visit(visited[i]);
    }
}

// breadth-first traversal starting from startLabel
void Graph::bfs(const string &startLabel, void visit(const string &label)) 
{
    if (!contains(startLabel)) 
    {
      return;
    }

    vector<string> visited;
    queue<Vertex *> queue;
    queue.push(vertexList[index("", startLabel)]);
    visited.emplace_back(queue.front()->value);

    while (!queue.empty()) 
    {
        Vertex *current = queue.front();
        queue.pop();

        for (int i = 0; i < current->edgeList.size(); i++) 
        {
            Vertex *neighbor = current->edgeList[i].first;
            if (find(visited.begin(), visited.end(), neighbor->value) == visited.end()) 
            {
                visited.emplace_back(neighbor->value);
                queue.push(neighbor);
            }
        }
    }

    string val;
    for (int i = 0; i < visited.size(); i++) 
    {
        val.append(visited[i]);
    }
    visit(val);
}

// store the weights in a map
// store the previous label in a map
std::pair<std::map<std::string, int>, std::map<std::string, std::string>>
Graph::dijkstra(const std::string &startLabel) const 
{
    
    const int tempLarge = 999999; 
    std::map<std::string, int> distanceMap;
    std::map<std::string, std::string> predecessorMap;

    if (contains(startLabel)) 
    {
        std::vector<std::string> visitedVertices;
        std::map<Vertex *, std::pair<std::vector<std::string>, int>> shortestPaths;
        std::map<Vertex *, int> tempDistances;

        distancePathHelper(vertexList[index("", startLabel)], shortestPaths, 
          0, visitedVertices, tempDistances, 0);
        shortestPaths.erase(vertexList.at(index("", startLabel)));

        while (shortestPaths.empty() == false) 
        {
            int minDistance = tempLarge; 
            int minIndex = 0;
            for (int i = 0; i < vertexList.size(); i++) {
                if (shortestPaths.count(vertexList[i]) != 0 &&
                    shortestPaths.at(vertexList[i]).second < minDistance) 
                {
                    minDistance = shortestPaths.at(vertexList[i]).second;
                    minIndex = i;
                }
            }
            distanceMap.emplace(vertexList[minIndex]->value, minDistance);
            predecessorMap.emplace(
                vertexList[minIndex]->value,
                shortestPaths.at(vertexList[minIndex])
                    .first.at(shortestPaths.at(vertexList[minIndex]).first.size() - 2));
            shortestPaths.erase(vertexList[minIndex]);
        }
        return std::make_pair(distanceMap, predecessorMap);
    }
    return std::make_pair(distanceMap, predecessorMap);
}


//////////////////////////////////////////////////////// extra credit algorithms ///////////////////////////////////////////////////////////

// extra credit: kruskal algorithm minimum spanning tree
// @return length of the minimum spanning tree or -1 if start vertex not

int Graph::mstKruskal(const string &startLabel,
    void visit(const string &from, const string &to, int weight)) const 
{
  
  vector<Vertex *> visited;
  int totalWeight = 0;
  int iterations = 0;
  // edge no greater than 
  const int tempLarge = 999999; 

  while (visited.size() < vertexList.size() && iterations < 100) 
  {
    int smallestEdgeWeight = tempLarge;
    int selectedVertexIndex = -1;
    int selectedEdgeIndex = -1;

    for (int i = 0; i < vertexList.size(); i++) 
    {
      for (int j = 0; j < vertexList[i]->edgeList.size(); j++) 
      {
        if (vertexList[i]->edgeList[j].second < smallestEdgeWeight &&
            find(visited.begin(), visited.end(), 
            vertexList[i]->edgeList[j].first) == visited.end()) 
        {
          iterations = 0;
          smallestEdgeWeight = vertexList[i]->edgeList[j].second;
          selectedVertexIndex = i;
          selectedEdgeIndex = j;
        }
      }
    }

    if (selectedVertexIndex != -1 && selectedEdgeIndex != -1) 
    {
      if (find(visited.begin(), visited.end(), 
      vertexList[selectedVertexIndex]->edgeList[selectedEdgeIndex].first) == visited.end()) 
      {
        visit(vertexList[selectedVertexIndex]->value,
        vertexList[selectedVertexIndex]->edgeList[selectedEdgeIndex].first->value,
        smallestEdgeWeight);
        visited.push_back(vertexList[selectedVertexIndex]->edgeList[selectedEdgeIndex].first);
        totalWeight += smallestEdgeWeight;
        if (find(visited.begin(), visited.end(), vertexList[selectedVertexIndex]) == visited.end()) 
        {
          visited.push_back(vertexList[selectedVertexIndex]);
        }
      }
    }
    iterations++;
  }
  return totalWeight;
}

// minimum spanning tree using Prim's algorithm. dijkstraPrimHelper gets
// the lowest path to every single vertex.
int Graph::mstPrim(const string &startLabel, 
void visit(const string &from, const string &to, int weight)) const 
{
    if (!contains(startLabel)) 
    {
        return -1; 
    }

    int totalWeight = 0;
    const int tempLarge = 999999; 
    map<Vertex*, pair<vector<string>, int>> shortestPaths;
    map<Vertex*, int> vertexDistances;

    distancePathHelper(vertexList[index("", startLabel)], shortestPaths, 0, {}, vertexDistances, 0);

    while (!shortestPaths.empty()) 
    {
        int smallestDistance = tempLarge;
        Vertex* smallestVertex = nullptr;
        for (const auto &entry : shortestPaths) 
        {
            if (entry.second.second < smallestDistance) 
            {
                smallestDistance = entry.second.second;
                smallestVertex = entry.first;
            }
        }

        if (smallestVertex != nullptr) 
        {
            const auto &path = shortestPaths[smallestVertex].first;
            if (path.size() > 1) 
            {
                visit(path[path.size() - 2], path.back(), vertexDistances[smallestVertex]);
                totalWeight += vertexDistances[smallestVertex];
            }
            shortestPaths.erase(smallestVertex);
        }
    }
    return totalWeight;
}

///////////////////////////////////////////////////////////// helper functions ////////////////////////////////////////////////////////////

void Graph::traverseVerticesHelper(const Vertex *curr, vector<string> &visited) 
{
    visited.push_back(curr->value); 
    for (int i = 0; i < curr->edgeList.size(); i++) 
    {
        const pair<Vertex*, int>& edge = curr->edgeList[i];
        if (find(visited.begin(), visited.end(), edge.first->value) == visited.end()) 
        {
            traverseVerticesHelper(edge.first, visited); 
        }
    }
}

void Graph::distancePathHelper(Vertex *currentVertex, map<Vertex *, 
    pair<vector<string>, int>> &shortestPaths, int currentDistance, 
    vector<string> currentPath, map<Vertex *, 
    int> &vertexDistances, int previousDistance) const
{
    if (shortestPaths.count(currentVertex) == 0 || 
    shortestPaths[currentVertex].second > currentDistance)
    {
        currentPath.push_back(currentVertex->value);
        shortestPaths[currentVertex] = make_pair(currentPath, currentDistance);
        vertexDistances[currentVertex] = previousDistance;
        
        for (vector<pair<Vertex *, int>>::const_iterator temp = 
        currentVertex->edgeList.begin(); 
        temp != currentVertex->edgeList.end(); temp++)
        {
            int newDistance = currentDistance + temp->second;
            distancePathHelper(temp->first, shortestPaths, 
            newDistance, currentPath, vertexDistances, temp->second);
        }
    }
}


