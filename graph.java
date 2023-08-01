ANS 1
import java.util.*;
 
// A class to store a graph edge
class Edge
{
    int source, dest;
 
    public Edge(int source, int dest)
    {
        this.source = source;
        this.dest = dest;
    }
}
 
// A class to represent a graph object
class Graph
{
    // A list of lists to represent an adjacency list
    List<List<Integer>> adjList = null;
 
    // Constructor
    Graph(List<Edge> edges, int n)
    {
        adjList = new ArrayList<>(n);
 
        for (int i = 0; i < n; i++) {
            adjList.add(i, new ArrayList<>());
        }
 
        // add edges to the undirected graph
        for (Edge edge: edges)
        {
            int src = edge.source;
            int dest = edge.dest;
 
            adjList.get(src).add(dest);
            adjList.get(dest).add(src);
        }
    }
}
 
class Main
{
 
    public static boolean DFS(Graph graph, int v, boolean[] discovered, int parent)
    {
     
        discovered[v] = true;
 
     
        for (int w: graph.adjList.get(v))
        {
           
            if (!discovered[w])
            {
                if (DFS(graph, w, discovered, v)) {
                    return true;
                }
            }
 
           
            else if (w != parent)
            {

		 return true;
            }
        }
 
        // No back-edges were found in the graph
        return false;
    }
 
    public static void main(String[] args)
    {
        // List of graph edges
        List<Edge> edges = Arrays.asList(
                        new Edge(0, 1), new Edge(0, 6), new Edge(0, 7),
                        new Edge(1, 2), new Edge(1, 5), new Edge(2, 3),
                        new Edge(2, 4), new Edge(7, 8), new Edge(7, 11),
                        new Edge(8, 9), new Edge(8, 10), new Edge(10, 11)
                        // edge (10, 11) introduces a cycle in the graph
                    );
 
       
        int n = 12;
 
    
        Graph graph = new Graph(edges, n);
 
      
        boolean[] discovered = new boolean[n];
 
     
        if (DFS(graph, 0, discovered, -1)) {
            System.out.println("The graph contains a cycle");
        }
        else {
            System.out.println("The graph doesn't contain any cycle");
        }
    }
}

ANS 2
	import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
 
// A class to store a graph edge
class Edge
{
    int source, dest;
 
    public Edge(int source, int dest)
    {
        this.source = source;
        this.dest = dest;
    }
}
 
// A class to represent a graph object
class Graph
{
    // A list of lists to represent an adjacency list
    List<List<Integer>> adjList = null;
 
    // Constructor
    Graph(List<Edge> edges, int n)
    {
        adjList = new ArrayList<>();
 
        for (int i = 0; i < n; i++) {
            adjList.add(new ArrayList<>());
        }
 
        // add edges to the directed graph
        for (Edge edge: edges) {
            adjList.get(edge.source).add(edge.dest);
        }
    }
}
 
class Main
{
    // Perform DFS on the graph and set the departure time of all
    // vertices of the graph
    private static int DFS(Graph graph, int v, boolean[] discovered,
                           int[] departure, int time)
    {
        // mark the current node as discovered
        discovered[v] = true;
 
        // do for every edge (v, u)
        for (int u: graph.adjList.get(v))
        {
            // if `u` is not yet discovered
            if (!discovered[u]) {
                time = DFS(graph, u, discovered, departure, time);
            }
        }
 
        // ready to backtrack
        // set departure time of vertex `v`
        departure[v] = time++;
 
        return time;
    }
 
    // Returns true if given directed graph is DAG
    public static boolean isDAG(Graph graph, int n)
    {
        // keep track of whether a vertex is discovered or not
        boolean[] discovered = new boolean[n];
 
        // keep track of the departure time of a vertex in DFS
        int[] departure = new int[n];
 
        int time = 0;
 
        // Perform DFS traversal from all undiscovered vertices
        // to visit all connected components of a graph
        for (int i = 0; i < n; i++)
        {
            if (!discovered[i]) {
                time = DFS(graph, i, discovered, departure, time);
            }
        }
 
        // check if the given directed graph is DAG or not
        for (int u = 0; u < n; u++)
        {
            // check if (u, v) forms a back-edge.
            for (int v: graph.adjList.get(u))
            {
                // If the departure time of vertex `v` is greater than equal
                // to the departure time of `u`, they form a back edge.
 
                // Note that departure[u] will be equal to
                // departure[v] only if `u = v`, i.e., vertex
                // contain an edge to itself
                if (departure[u] <= departure[v]) {
                    return false;
                }
            }
        }
 
        // no back edges
        return true;
    }
 
    public static void main(String[] args)
    {
        // List of graph edges as per the above diagram
        List<Edge> edges = Arrays.asList(
                new Edge(0, 1), new Edge(0, 3), new Edge(1, 2),
                new Edge(1, 3), new Edge(3, 2), new Edge(3, 4),
                new Edge(3, 0), new Edge(5, 6), new Edge(6, 3)
        );
 
        // total number of nodes in the graph (labelled from 0 to 6)
        int n = 7;
 
        // build a graph from the given edges
        Graph graph = new Graph(edges, n);
 
        // check if the given directed graph is DAG or not
        if (isDAG(graph, n)) {
            System.out.println("The graph is a DAG");
        }
        else {
            System.out.println("The graph is not a DAG");
        }
    }
}

ANS 3
import java.util.*;
 
// A class to store a graph edge
class Edge
{
    public final int src, dest, weight;
 
    private Edge(int src, int dest, int weight)
    {
        this.src = src;
        this.dest = dest;
        this.weight = weight;
    }
 
    // Factory method for creating an immutable instance of `Edge`
    public static Edge of(int a, int b, int c) {
        return new Edge(a, b, c);        // calls private constructor
    }
}
 
// A BFS Node
class Node
{
    int vertex, depth, weight;
 
    Node(int vertex, int depth, int weight)
    {
        this.vertex = vertex;
        this.depth = depth;
        this.weight = weight;
    }
}
 
// A class to represent a graph object
class Graph
{
    // A list of lists to represent an adjacency list
    List<List<Edge>> adjList = new ArrayList<>();
 
    // Graph Constructor
    public Graph(List<Edge> edges, int n)
    {
        // resize the list to `n` elements of type `List<Edge>`
        for (int i = 0; i < n; i++) {
            adjList.add(new ArrayList<>());
        }
 
        // add edges to the directed graph
        for (Edge e: edges) {
            adjList.get(e.src).add(e);
        }
    }
}
 
class Main
{
    // Perform BFS on graph `g` starting from vertex `v`
    public static int findLeastCost(Graph g, int src, int dest, int m)
    {
        // create a queue for doing BFS
        Queue<Node> q = new ArrayDeque<>();
 
        // enqueue source vertex
        q.add(new Node(src, 0, 0));
 
        // stores least-cost from source to destination
        int minCost = Integer.MAX_VALUE;
 
        // loop till queue is empty
        while (!q.isEmpty())
        {
            // dequeue front node
            Node node = q.poll();
 
            int v = node.vertex;
            int depth = node.depth;
            int cost = node.weight;
 
            // if the destination is reached and BFS depth is equal to `m`,
            // update the minimum cost calculated so far
            if (v == dest && depth == m) {
                minCost = Math.min(minCost, cost);
            }
 
            // don't consider nodes having a BFS depth more than `m`.
            // This check will result in optimized code and handle cycles
            // in the graph (otherwise, the loop will never break)
            if (depth > m) {
                break;
            }
 
            // do for every adjacent edge of `v`
            for (Edge edge: g.adjList.get(v))
            {
                // push every vertex (discovered or undiscovered) into
                // the queue with depth as +1 of parent and cost equal
                // to the cost of parent plus the current edge weight
                q.add(new Node(edge.dest, depth + 1, cost + edge.weight));
            }
        }
 
        // return least-cost
        return minCost;
    }
 
    public static void main(String[] args)
    {
        // List of graph edges as per the above diagram
        List<Edge> edges = Arrays.asList(
                        Edge.of(0, 6, -1), Edge.of(0, 1, 5), Edge.of(1, 6, 3),
                        Edge.of(1, 5, 5), Edge.of(1, 2, 7), Edge.of(2, 3, 8),
                        Edge.of(3, 4, 10), Edge.of(5, 2, -1), Edge.of(5, 3, 9),
                        Edge.of(5, 4, 1), Edge.of(6, 5, 2), Edge.of(7, 6, 9),
                        Edge.of(7, 1, 6));
 
        // total number of nodes in the graph (labelled from 0 to 7)
        int n = 8;
 
        // build a graph from the given edges
        Graph g = new Graph(edges, n);
 
        int src = 0, dest = 3, m = 4;
 
        // Perform modified BFS traversal from source vertex `src`
        System.out.print(findLeastCost(g, src, dest, m));
    }
}

ANS 4

	import java.util.*;
 
// A class to store a graph edge
class Edge
{
    int source, dest;
 
    public Edge(int source, int dest)
    {
        this.source = source;
        this.dest = dest;
    }
}
 
// A BFS Node
class Node
{
    // stores current vertex number and the current depth of
    // BFS (how far away from source current node is)
    int vertex, depth;
 
    public Node(int vertex, int depth)
    {
        this.vertex = vertex;
        this.depth = depth;
    }
}
 
// A class to represent a graph object
class Graph
{
    // A list of lists to represent an adjacency list
    List<List<Integer>> adjList = null;
 
    // Constructor
    Graph(List<Edge> edges, int n)
    {
        adjList = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            adjList.add(new ArrayList<>());
        }
 
        // add edges to the directed graph
        for (Edge edge: edges) {
            adjList.get(edge.source).add(edge.dest);
        }
    }
}
 
class Main
{
    // Perform BFS on graph `graph` starting from vertex `v`
    public static int findTotalPaths(Graph graph, int src, int dest, int m)
    {
        // create a queue for doing BFS
        Queue<Node> q = new ArrayDeque<>();
 
        // enqueue source vertex
        q.add(new Node(src, 0));
 
        // stores number of paths from source to destination
        // having exactly `m` edges
        int count = 0;
 
        // loop till queue is empty
        while (!q.isEmpty())
        {
            // dequeue front node
            Node node = q.poll();
 
            int v = node.vertex;
            int depth = node.depth;
 
            // if the destination is reached and BFS depth is equal to `m`,
            // update count
            if (v == dest && depth == m) {
                count++;
            }
 
            // don't consider nodes having a BFS depth more than `m`.
            // This check will result in optimized code and handle cycles
            // in the graph (otherwise, the loop will never break)
            if (depth > m) {
                break;
            }
 
            // do for every adjacent vertex `u` of `v`
            for (int u: graph.adjList.get(v))
            {
                // push every vertex (discovered or undiscovered) into the queue
                q.add(new Node(u, depth + 1));
            }
        }
 
        // return number of paths from source to destination
        return count;
    }
 
    public static void main(String[] args)
    {
        // List of graph edges as per the above diagram
        List<Edge> edges = Arrays.asList(
                new Edge(0, 6), new Edge(0, 1), new Edge(1, 6), new Edge(1, 5),
                new Edge(1, 2), new Edge(2, 3), new Edge(3, 4), new Edge(5, 2),
                new Edge(5, 3), new Edge(5, 4), new Edge(6, 5), new Edge(7, 6),
                new Edge(7, 1));
 
        // total number of nodes in the graph
        int n = 8;
 
        // construct graph
        Graph graph = new Graph(edges, n);
 
        int src = 0, dest = 3, m = 4;
 
        // Do modified BFS traversal from the source vertex src
        System.out.println(findTotalPaths(graph, src, dest, m));
    }
}

ANS 5import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
 
// A class to store a graph edge
class Edge
{
    int source, dest;
 
    public Edge(int source, int dest)
    {
        this.source = source;
        this.dest = dest;
    }
}
 
// A class to represent a graph object
class Graph
{
    // A list of lists to represent an adjacency list
    List<List<Integer>> adjList;
 
    // stores indegree of a vertex
    List<Integer> in;
 
    // Constructor
    Graph(int n)
    {
        adjList = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            adjList.add(new ArrayList<>());
        }
 
        // initialize indegree of each vertex by 0
        in = new ArrayList<>(Collections.nCopies(n, 0));
    }
 
    // Utility function to add an edge (u, v) to the graph
    void addEdge(int u, int v)
    {
        // add an edge from source to destination
        adjList.get(u).add(v);
 
        // increment in-degree of destination vertex by 1
        in.set(v, in.get(v) + 1);
    }
}
 
class Main
{
    // Function to perform DFS traversal on the graph
    public static void DFS(Graph graph, int v, boolean[] discovered)
    {
        // mark the current node as discovered
        discovered[v] = true;
 
        // do for every edge (v, u)
        for (int u: graph.adjList.get(v))
        {
            // `u` is not discovered
            if (!discovered[u]) {
                DFS(graph, u, discovered);
            }
        }
    }
 
    // Function to replace all the directed edges of the graph with undirected edges
    // and produce an undirected graph
    public static Graph getUndirectedGraph(Graph graph, int n)
    {
        Graph g = new Graph(n);
        for (int u = 0; u < n; u++)
        {
            for (int v: graph.adjList.get(u)) {        // for every edge (u, v)
                g.addEdge(v, u);
                g.addEdge(u, v);
            }
        }
        return g;
    }
 
    // Function to check if all vertices with a non-zero degree in a graph
    // belong to a single connected component
    public static boolean isConnected(Graph graph, int n)
    {
        // keep track of all previously visited vertices in DFS
        boolean[] visited = new boolean[n];
 
        // start DFS from the first vertex with a non-zero degree
        for (int i = 0; i < n; i++)
        {
            if (graph.adjList.get(i).size() > 0)
            {
                DFS(graph, i, visited);
                break;
            }
        }
 
        // if a single DFS call couldn't visit all vertices with a non-zero degree,
        // the graph contains more than one connected component
        for (int i = 0; i < n; i++)
        {
            if (!visited[i] && graph.adjList.get(i).size() > 0) {
                return false;
            }
        }
 
        return true;
    }
 
    // Function to check if a directed graph has an Eulerian path
    public static boolean hasEulerPath(Graph graph, int n)
    {
        /*
            The following loop checks the following conditions to determine if an
            Eulerian path can exist or not:
                a. At most one vertex in the graph has `out-degree = 1 + in-degree`.
                b. At most one vertex in the graph has `in-degree = 1 + out-degree`.
                c. Rest all vertices have `in-degree == out-degree`.
 
            If either of the above condition fails, the Euler path can't exist.
        */
 
        boolean x = false, y = false;
        for (int i = 0; i < n; i++)
        {
            int out_degree = graph.adjList.get(i).size();
            int in_degree = graph.in.get(i);
 
            if (out_degree != in_degree)
            {
                if (!x && out_degree - in_degree == 1) {
                    x = true;
                }
                else if (!y && in_degree - out_degree == 1) {
                    y = true;
                }
                else {
                    return false;
                }
            }
        }
 
        // consider given edges as undirected and check if all non-isolated vertices
        // belong to a single connected component
        return isConnected(getUndirectedGraph(graph, n), n);
    }
 
    public static void main(String[] args)
    {
        // List of graph edges as per the above diagram
        List<Edge> edges = Arrays.asList(new Edge(0, 1), new Edge(1, 2),
                new Edge(2, 3), new Edge(3, 1), new Edge(1, 4), new Edge(4, 3),
                new Edge(3, 0), new Edge(0, 5), new Edge(5, 4)
        );
 
        // total number of nodes in the graph (labelled from 0 to 5)
        int n = 6;
 
        // build a directed graph from the above edges
        Graph graph = new Graph(n);
 
        // add edges to the directed graph
        for (Edge edge: edges) {
            graph.addEdge(edge.source, edge.dest);
        }
 
        if (hasEulerPath(graph, n)) {
            System.out.println("The graph has an Eulerian path");
        }
        else {
            System.out.println("The Graph doesn't have an Eulerian path");
        }
    }
}



									ASSIGNMENT 2
ANS 1
	import java.util.*;
 
// A class to store a graph edge
class Edge
{
    int source, dest;
 
    public Edge(int source, int dest)
    {
        this.source = source;
        this.dest = dest;
    }
}
 
// A class to represent a graph object
class Graph
{
    // A list of lists to represent an adjacency list
    List<List<Integer>> adjList = null;
 
    // Total number of nodes in the graph
    int n;
 
    // Constructor
    Graph(List<Edge> edges, int n)
    {
        this.adjList = new ArrayList<>();
        this.n = n;
 
        for (int i = 0; i < n; i++) {
            adjList.add(new ArrayList<>());
        }
 
        // add edges to the undirected graph
        for (Edge edge: edges)
        {
            int src = edge.source;
            int dest = edge.dest;
 
            // add an edge from source to destination
            adjList.get(src).add(dest);
 
            // add an edge from destination to source
            adjList.get(dest).add(src);
        }
    }
}
 
class Main
{
    // Perform BFS on the graph starting from vertex `v`
    public static boolean isBipartite(Graph graph)
    {
        // get total number of nodes in the graph
        int n = graph.n;
 
        // start from any node as the graph is connected and undirected
        int v = 0;
 
        // to keep track of whether a vertex is discovered or not
        boolean[] discovered = new boolean[n];
 
        // stores the level of each vertex in BFS
        int[] level = new int[n];
 
        // mark the source vertex as discovered and
        // set its level to 0
        discovered[v] = true;
        level[v] = 0;
 
        // create a queue to do BFS and enqueue
        // source vertex in it
        Queue<Integer> q = new ArrayDeque<>();
        q.add(v);
 
        // loop till queue is empty
        while (!q.isEmpty())
        {
            // dequeue front node and print it
            v = q.poll();
 
            // do for every edge (v, u)
            for (int u: graph.adjList.get(v))
            {
                // if vertex `u` is explored for the first time
                if (!discovered[u])
                {
                    // mark it as discovered
                    discovered[u] = true;
 
                    // set level one more than the level of the parent node
                    level[u] = level[v] + 1;
 
                    // enqueue vertex
                    q.add(u);
                }
                // if the vertex has already been discovered and the
                // level of vertex `u` and `v` are the same, then the
                // graph contains an odd-cycle and is not bipartite
                else if (level[v] == level[u]) {
                    return false;
                }
            }
        }
 
        return true;
    }
 
    public static void main(String[] args)
    {
        // List of graph edges
        List<Edge> edges = Arrays.asList(
                    new Edge(0, 1), new Edge(1, 2), new Edge(1, 7), new Edge(2, 3),
                    new Edge(3, 5), new Edge(4, 6), new Edge(4, 8), new Edge(7, 8)
                    // if we add edge (1, 3), the graph becomes non-bipartite
                );
 
        // total number of nodes in the graph (0 to 8)
        int n = 9;
 
        // build a graph from the given edges
        Graph graph = new Graph(edges, n);
 
        if (isBipartite(graph)) {
            System.out.println("Graph is bipartite");
        }
        else {
            System.out.println("Graph is not bipartite");
        }
    }
}

ANS 2
		import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;
 
class Main
{
    // Function to print the itinerary starting from a given source `src`
    private static void printItinerary(Map<String, String> map, String src, int flag)
    {
        String dest = map.get(src);
        if (dest == null) {
            return;
        }
      
        if(flag == 0)
            System.out.print(src + " —> " + dest);
        else 
            System.out.print(" —> " + dest);
        printItinerary(map, dest, 1);
    }
 
    // Function to find the itinerary from the given list of departure
    // and arrival airports
    private static void findItinerary(String[][] input)
    {
        // construct a map from the given list of tickets (departure —> arrival)
        Map<String, String> tickets = Stream.of(input)
                .collect(Collectors.toMap(p -> p[0], p -> p[1]));
 
        // construct a set of destination airports
        Set<String> destinations = new HashSet<>(tickets.values());
 
        // consider each departure airport to find the source airport
        for (String airport: tickets.keySet())
        {
            // source airport will not be present in the list of destination airports
            if (!destinations.contains(airport))
            {
                // when the source airport is found, print the itinerary
                printItinerary(tickets, airport, 0);
                return;
            }
        }
    }
 
    public static void main(String[] args)
    {
        // input: list of tickets
        String[][] input = new String[][]{
                {"LAX", "DXB"},
                {"DFW", "JFK"},
                {"LHR", "DFW"},
                {"JFK", "LAX"}
        };
 
        findItinerary(input);
    }
}

ANS 3
	import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;
 
class Main
{
    // Function to print the itinerary starting from a given source `src`
    private static void printItinerary(Map<String, String> map, String src, int flag)
    {
        String dest = map.get(src);
        if (dest == null) {
            return;
        }
      
        if(flag == 0)
            System.out.print(src + " —> " + dest);
        else 
            System.out.print(" —> " + dest);
        printItinerary(map, dest, 1);
    }
 
    // Function to find the itinerary from the given list of departure
    // and arrival airports
    private static void findItinerary(String[][] input)
    {
        // construct a map from the given list of tickets (departure —> arrival)
        Map<String, String> tickets = Stream.of(input)
                .collect(Collectors.toMap(p -> p[0], p -> p[1]));
 
        // construct a set of destination airports
        Set<String> destinations = new HashSet<>(tickets.values());
 
        // consider each departure airport to find the source airport
        for (String airport: tickets.keySet())
        {
            // source airport will not be present in the list of destination airports
            if (!destinations.contains(airport))
            {
                // when the source airport is found, print the itinerary
                printItinerary(tickets, airport, 0);
                return;
            }
        }
    }
 
    public static void main(String[] args)
    {
        // input: list of tickets
        String[][] input = new String[][]{
                {"LAX", "DXB"},
                {"DFW", "JFK"},
                {"LHR", "DFW"},
                {"JFK", "LAX"}
        };
 
        findItinerary(input);
    }
}

ANS 4
	import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
 
// A class to store a graph edge
class Edge
{
    int source, dest;
 
    public Edge(int source, int dest)
    {
        this.source = source;
        this.dest = dest;
    }
}
 
// A class to represent a graph object
class Graph
{
    // A list of lists to represent an adjacency list
    List<List<Integer>> adjList = null;
 
    // Constructor
    Graph(List<Edge> edges, int n)
    {
        adjList = new ArrayList<>();
 
        for (int i = 0; i < n; i++) {
            adjList.add(new ArrayList<>());
        }
 
        // add edges to the directed graph
        for (Edge edge: edges) {
            adjList.get(edge.source).add(edge.dest);
        }
    }
}
 
class Main
{
    // Function to perform DFS traversal on the graph on a graph
    public static int DFS(Graph graph, int v, boolean[] discovered,
                    int[] arrival, int[] departure, int time)
    {
        // set the arrival time of vertex `v`
        arrival[v] = ++time;
 
        // mark vertex as discovered
        discovered[v] = true;
 
        for (int i: graph.adjList.get(v))
        {
            if (!discovered[i]) {
                time = DFS(graph, i, discovered, arrival, departure, time);
            }
        }
 
        // set departure time of vertex `v`
        departure[v] = ++time;
 
        return time;
    }
 
    public static void main(String[] args)
    {
        // List of graph edges as per the above diagram
        List<Edge> edges = Arrays.asList(
                new Edge(0, 1), new Edge(0, 2), new Edge(2, 3), new Edge(2, 4),
                new Edge(3, 1), new Edge(3, 5), new Edge(4, 5), new Edge(6, 7)
        );
 
        // total number of nodes in the graph (labelled from 0 to 7)
        int n = 8;
 
        // build a graph from the given edges
        Graph graph = new Graph(edges, n);
 
        // array to store the arrival time of vertex
        int[] arrival = new int[n];
 
        // array to store the departure time of vertex
        int[] departure = new int[n];
 
        // mark all the vertices as not discovered
        boolean[] discovered = new boolean[n];
        int time = -1;
 
        // Perform DFS traversal from all undiscovered nodes to
        // cover all unconnected components of a graph
        for (int i = 0; i < n; i++)
        {
            if (!discovered[i]) {
                time = DFS(graph, i, discovered, arrival, departure, time);
            }
        }
 
        // print arrival and departure time of each vertex in DFS
        for (int i = 0; i < n; i++)
        {
            System.out.println("Vertex " + i + " (" + arrival[i] + ", " +
                departure[i] + ")");
        }
    }
}

ANS 5
		import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
 
// A class to store a graph edge
class Edge
{
    int source, dest;
 
    public Edge(int source, int dest)
    {
        this.source = source;
        this.dest = dest;
    }
}
 
class Graph
{
    List<List<Integer>> adjList;
 
    Graph(List<Edge> edges, int n)
    {
        adjList = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            adjList.add(new ArrayList<>());
        }
 
        // add edges to the directed graph
        for (Edge edge: edges) {
            adjList.get(edge.source).add(edge.dest);
        }
    }
}
 
class Main
{
   
    public static void DFS(Graph graph, int v, boolean[] discovered)
    {
        
        discovered[v] = true;
 
       
        for (int u: graph.adjList.get(v))
        {
           
            if (!discovered[u]) {
                DFS(graph, u, discovered);
            }
        }
    }
 
   
    public static int findRootVertex(Graph graph, int n)
    {
        // to keep track of all previously visited vertices in DFS
        boolean[] visited = new boolean[n];
 
        // find the last starting vertex `v` in DFS
        int v = 0;
        for (int i = 0; i < n; i++)
        {
            if (!visited[i])
            {
                DFS(graph, i, visited);
                v = i;
            }
        }
 
        
        Arrays.fill(visited, false);
 
       
        DFS(graph, v, visited);
 
       
        for (int i = 0; i < n; i++)
        {
            if (!visited[i]) {
                return -1;
            }
        }
 
        
        return v;
    }
 
    public static void main(String[] args)
    {
      
        List<Edge> edges = Arrays.asList(new Edge(0, 1), new Edge(1, 2),
                    new Edge(2, 3), new Edge(3, 0), new Edge(4, 3),
                    new Edge(4, 5), new Edge(5, 0));
 
       
        int n = 6;
 
       
        Graph graph = new Graph(edges, n);
 
       
        int root = findRootVertex(graph, n);
 
        if (root != -1) {
            System.out.println("The root vertex is " + root);
        }
        else {
            System.out.println("The root vertex does not exist");
        }
    }
}