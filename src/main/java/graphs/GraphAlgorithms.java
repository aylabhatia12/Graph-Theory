/**
 * 
 * 
 * NAME: Ayla Bhatia
 */

package graphs;

import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;

import java.util.Deque; // addFirst = push, removeFirst = pop, getFirst = peek
import java.util.Collections;

/** 
 * Suite of graph-based algorithms. 
 */
public class GraphAlgorithms {

  /**
   * Performs a recursive depth-first traversal of a directed graph
   * starting at a source vertex.
   * @param g The graph to search
   * @param s The source vertex to search from
   * @param digraph If false, treat g as undirected
   * @returns A search tree as a map with vertices as keys and
   * search-tree parents as corresponding values, or null if the
   * source vertex is invalid.
   */
  public static Map<Integer,Integer> dfs_rec(Graph g, int s, boolean digraph) {
    // TODO: implement this function
    //Returning null if the starting vertex does not exist
    if (!g.hasVertex(s)) {
        return null; 
    }

    boolean[] visitedArray = new boolean[g.vertices()];   // A boolean array to store if each vertex has been visited
    Map<Integer, Integer> parentMap = new HashMap<>(); // Map that will store the parents
   
   parentMap.put(s, -1); // Setting the parent of the starting node
    // Start DFS by calling the recursive function dfs_visit
    dfs_visit(g, s, visitedArray, parentMap, digraph);
    return parentMap;
  }

  // Helper method for the recursive part of dfs_rec, that calls itself until all vertices are visited
  private static Map<Integer,Integer> dfs_visit(Graph g, int current, boolean[] visited,
                               Map<Integer, Integer> parents , boolean digraph) {
        // Mark the current node as visited
    visited[current] = true; 
    // The parent of the current node is the key in parentMap, 
    // we will check it when we start visiting neighbors

    // Get neighbors into a set based on whether the graph is directed or undirected
    Set<Integer> neighbors;
    if (digraph) {
        neighbors = g.out(current); 
    } else {
        neighbors = g.adj(current); 
    }

    // Look at all neighbors    
    for (int neighbor : neighbors) {
        if (!visited[neighbor]) {
            // Set the parent of the neighbor to be the current vertex
            parents.put(neighbor, current);
            // Recursive call to visit the neighbor
            dfs_visit(g, neighbor, visited, parents, digraph); 
        } 
    }    
    return parents;
}
  
  /**
   * Performs an interative depth-first traversal of a directed graph
   * starting at a source vertex.
   * @param g The graph to search
   * @param s The source vertex to search from
   * @param digraph If false, treat g as undirected
   * @returns A search tree as a map with vertices as keys and
   * search-tree parents as corresponding values, or null if the
   * source vertex is invalid.
   */
  public static Map<Integer,Integer> dfs(Graph g, int s, boolean digraph) {
    // TODO: implement this function
    if (!g.hasVertex(s)) {
          return null;
      }

      // Map to store parents
      Map<Integer, Integer> parentMap = new HashMap<>();
      // USing a set to track visited vertices
      Set<Integer> visited = new HashSet<>();
      // Stack for DFS
      Deque<Integer> stack = new LinkedList<>();

      // Add source vertex with no parent (-1)
      parentMap.put(s, -1);
      stack.addFirst(s);
      
      while (!stack.isEmpty()) {
          int current = stack.removeFirst();// Get the current vertex from the stack

          if (visited.contains(current)) { // If current vertext has already been visited, do nothing and continue looking at the stack.
              continue;
          }

          visited.add(current); // If current veterx has not been visited, add it to the visited set now.

          // Get neighbors of the current vertex into a set based on whether graph is directed or undirected
          Set<Integer> neighbors;
          if (digraph) {
              neighbors = g.out(current); 
          } else {
              neighbors = g.adj(current); 
          }

          // Sort neighbors in descending order to match recursive dfs. 
          // I had to do this after comparing the output from the recursive dfs method and getting different result in the iterative dfs.
          List<Integer> sortedNeighbors = new ArrayList<>(neighbors);
          Collections.sort(sortedNeighbors, Collections.reverseOrder());

          //Now look at all neighbors of the current vertex
          for (int neighbor : sortedNeighbors) {
              if (!visited.contains(neighbor)) { // If a neighbor has not been visited yet
                  // Add neighbor to stack. Set its paren, which is the current vertex, in dfs tree
                  stack.addFirst(neighbor);
                  parentMap.put(neighbor, current);
              }
          }
      }
      return parentMap;
  }
  /**
   * Checks if the given directed graph is acyclic using a modified,
   * iterative version of DFS. 
   * @param g The graph to check
   * @returns True if the graph does not have any cycles, and false
   * otherwise. 
   */  
 public static boolean acyclic(Graph g) {
    boolean[] visited = new boolean[g.vertices()]; // boolean array to track visited vertices
    boolean[] recStack = new boolean[g.vertices()]; // boolean array to track vertices in the recursion stack
    Deque<Integer> stack = new LinkedList<>(); // Stack for DFS

    // Perform DFS from each vertex
    for (int i = 0; i < g.vertices(); i++) {
        if (!visited[i]) { //If the vetex has not been visited yet
            // Start DFS from the current vertex
            stack.addFirst(i);
            while (!stack.isEmpty()) {
                int vertex = stack.getFirst(); // Get the vertex on top of the stack

                if (!visited[vertex]) { // If the vertes from the top of the stack has not been visited yet
                    visited[vertex] = true; // Mark as visited
                    recStack[vertex] = true; // Add to recursion stack

                    // Get its neighbors and check for cycles
                    List<Integer> neighbors = new ArrayList<>(g.out(vertex));
                    boolean allNeighborsChecked = true; // a variable to mark if all neighbors of the vertex have been checked

                    for (int neighbor : neighbors) {
                        if (!visited[neighbor]) {
                            stack.addFirst(neighbor); // Push unvisited neighbor to the dfs stack
                            allNeighborsChecked = false; // There are unvisited neighbors
                        } 
                        else if (recStack[neighbor]) {
                            return false; // It's a cycle if the neighbor is in the recursion stack
                        }
                    }

                    // If all neighbors are checked, remove the vertex from recursion stack
                    if (allNeighborsChecked) {
                        stack.removeFirst();
                        recStack[vertex] = false; // Remove from recursion stack
                    }
                } 
                else {
                    // If already visited, just remove from stack and remove from recursion stack
                    stack.removeFirst();
                    recStack[vertex] = false; // Remove from recursion stack
                }
            }
        }
    }
    return true; // No cycles detected
 }

  /**
   * Checks if the given undirected graph is bipartite by using a
   * modified, iterative version of DFS that tries to perform a
   * 2-coloring.
   * @param g The graph to check
   * @returns True if the graph is bipartite, and false
   * otherwise. 
   */
public static boolean bipartite(Graph g) {
    int[] color = new int[g.vertices()]; // Int arraye to store color for each vertext. 0-unvisited, 1-color1, 2-color2
    Deque<Integer> stack = new LinkedList<>();

    //Look at each vertex
    for (int i = 0; i < g.vertices(); i++) {
        if (color[i] == 0) { // if it is Unvisited
            color[i] = 1; // mark it color 1
            stack.addFirst(i); // Add to the stack

            while (!stack.isEmpty()) {
                int vertex = stack.removeFirst(); // Get the top vertex from the stack

                // Check all neighbors of the vertex
                for (int neighbor : g.adj(vertex)) {
                    if (color[neighbor] == 0) { // If a neighbor is unvisited
                        color[neighbor] = 3 - color[vertex]; // mark the neighbor with the alternate color
                        stack.addFirst(neighbor); // Push neighbor to stack
                    } 
                    else if (color[neighbor] == color[vertex]) { // If the color of the neighbor anfd color of the top vertex are the same, it's not bipartite graph.
                        return false; // Not Bipartitie
                    }
                }
            }
        }
    }
    return true; // Bipartite
}
  
  /**
   * Computes a topological sort of the given directed graph, if one
   * exists,  using a modified version of iterative DFS.
   * @param g The graph to sort
   * @returns A list giving a topological sort of the graph vertices,
   * if one exists, otherwise returns an empty list. 
   */
public static List<Integer> topologicalSort(Graph g) {
    // Result list to store the topological sorted order
      List<Integer> result = new ArrayList<>();

      // If graph is empty(no vertices), return empty list
      // Also check If graph has cycles, topological sorting is not possible, return emplty list
      if (g.vertices() == 0 || !acyclic(g)) {
          return result;
      }

      boolean[] visited = new boolean[g.vertices()]; //Array to track visited vertices(
      boolean[] onStack = new boolean[g.vertices()]; //Array to track vertices currently in the recursion stack for maintingin order

      Deque<Integer> stack = new LinkedList<>();// Stack for DFS
      Deque<Integer> orderStack = new LinkedList<>(); // Stack to store the topological ordering

      // Go through all vertices of the graph
      for (int vertex = 0; vertex < g.vertices(); vertex++) {
          if (!visited[vertex]) {
              stack.addFirst(vertex); // add to the dfs stack if the vertext has not been visited

              while (!stack.isEmpty()) {
                  int current = stack.getFirst();

                  if (!visited[current]) { // If current vertex is not visited
                      visited[current] = true; //mark it as visited
                      onStack[current] = true;  //put in on recursion stack
                  } 
                  else 
                        if (onStack[current]) { //if the vertex has been visited and is on the recursion stack, that shows all its children have been visited.
                            stack.removeFirst(); //remove from the dfs stack
                            onStack[current] = false;
                            orderStack.addFirst(current); // add it to the orderStack
                            continue; 
                        } 
                        else //If the vertex has been visited and is NOT on the recursion stack means already fully explored this vertex
                        {                            
                            stack.removeFirst();
                            continue;
                        }

                        // Add unvisited neighbors to stack
                        Set<Integer> neighbors = g.out(current);
                        boolean allVisited = true;

                        for (int neighbor : neighbors) {
                            if (!visited[neighbor]) {
                                stack.addFirst(neighbor);
                                allVisited = false;
                            } else if (onStack[neighbor]) {
                                // Skip vertices already in current path
                                continue;
                            }
                        }

                        // If all neighbors are visited, we're done with this vertex
                        if (allVisited) {
                            stack.removeFirst();
                            onStack[current] = false;
                            orderStack.addFirst(current);
                        }
              }
          }
      }
      // Transfer from orderStack to result list
      while (!orderStack.isEmpty()) {
          result.add(orderStack.removeFirst());
      }
      return result;
  }  

  /**
   * Performs a breadth-first traversal of a directed graph starting
   * at a source vertex.
   * @param g The graph to search
   * @param s The source vertex to search from
   * @param digraph If false, treat g as undirected
   * @returns A search tree as a map with vertices as keys and
   * search-tree parents as corresponding values, or null if the
   * source vertex is invalid. 
   */
    public static Map<Integer,Integer> bfs(Graph g, int s, boolean digraph) {
        // @TODO: implement this function
            if (g.hasVertex(s) == false) {
                return null; // source vertex does not exist
            }

            Map<Integer, Integer> searchTreeMap = new HashMap<>(); // To store search tree vertices
            boolean[] visitedVertices = new boolean[g.vertices()]; // Boolean array for visited vertices with the size equal to total number of vertices 
            Queue<Integer> queue = new LinkedList<>(); // Queue for BFS

            visitedVertices[s] = true; // Mark source vertex as visited
            queue.add(s); // Add the starting vertex to the bfs queue
            searchTreeMap.put(s, -1); // Add source vertex to the search tree with no parent

            while (queue.isEmpty() == false) { //If the bfs queue is not empty
                int current = queue.remove(); //remove the first vertex
                // Go through all outgoing neighbors
                Set<Integer> neighbors = g.out(current); // get all out edges for the current vertex
                for (Integer neighbor : neighbors) { // check all out edges
                    if (visitedVertices[neighbor] == false) { //This neighbor has not been visited
                        visitedVertices[neighbor] = true; // Mark this neighbor vertex as visited
                        queue.add(neighbor); // Add this neighbor vertex to the bfs queue
                        searchTreeMap.put(neighbor, current); // Add thus neighbor and its parent to the search map
                    }
                }

                //Chgecking for incoming edges if it's an undirected graph
                if (digraph == false) {
                    neighbors = g.in(current);
                    for (Integer neighbor : neighbors) {
                        if (visitedVertices[neighbor] == false) { //This neighbor has not been visited
                            visitedVertices[neighbor] = true; // Mark this neighbor vertex as visited
                            queue.add(neighbor); // Add this neighbor vertex to the bfs queue
                            searchTreeMap.put(neighbor, current); // Add thus neighbor and its parent to the search map
                        }
                    }
                }
            }
            return searchTreeMap; // Return the bfs search tree
        }

  /**
   * Checks if a source (src) vertex can reach a destination (dst)
   * vertex (i.e., dst is reachable from src) in a graph using a
   * modified version of the BFS algorithm. Assumes all vertices can
   * trivially reach themselves (even if no self edge).
   * @param g The graph to search
   * @param src The source vertex to search from
   * @param dst The destination vertex (for a path from src to dst)
   * @param digraph If false, treat g as undirected
   * @returns The shortest path from src to dst as a sequence of vertices
   *          (forming the path) or the empty sequence if dst is not
   *          reachable from src
   */
  public static List<Integer> reachable(Graph g, int src, int dst, boolean digraph) {
    // @TODO: implement this function
    // NOTE: You must not call your BFS function to implement reachable
    // Check if src and dest are valid
        if (g.hasVertex(src)  == false || g.hasVertex(dst) == false) {
            return new ArrayList<>(); // Return an empty list if src or dst vertices do not exist
        }

        List<Integer> path = new ArrayList<>(); // To store the path from src to dst
        Map<Integer, Integer> parentMap = new HashMap<>(); // To build the path
        boolean[] visitedVertices = new boolean[g.vertices()]; // Boolean array for visited vertices
        Queue<Integer> queue = new LinkedList<>(); // Queue for BFS

        visitedVertices[src] = true; // Mark source vertex as visited
        queue.add(src); // Add the starting vertex to the bfs queue
        parentMap.put(src, null); // Set the parent of the scource vertex as null, no parent

        while (queue.isEmpty() == false) { //Loop until the queue is empty
            int current = queue.poll(); // Remove the first vertex from the queue

            // stop if the vertex removed from the queue is the destinatiuon vertex
            if (current == dst) {
                break;
            }

            // Go through all outgoing neighbors
            Set<Integer> neighbors = g.out(current); //get all out edges for the current vertex
            for (Integer neighbor : neighbors) {
                if (!visitedVertices[neighbor]) { // if this neighbor has not been visited
                    visitedVertices[neighbor] = true; // Mark this neighbor veretx as visited now
                    queue.add(neighbor); // add this vertex to the queue
                    parentMap.put(neighbor, current); //add this vertes and parent to the map
                }
            }

            
            if (!digraph) { // If the graph is undirected
                neighbors = g.in(current); // get all in edges for the current vertex
                for (Integer neighbor : neighbors) {
                    if (!visitedVertices[neighbor]) { //if this neighbor has not been visited
                        visitedVertices[neighbor] = true; // Mark this neighbor veretx as visited now
                        queue.add(neighbor); //add this vertex to the queue
                        parentMap.put(neighbor, current); //add this vertes and parent to the map
                    }
                }
            }
        }

        // build the path from src to dst
        Integer current = dst;
        while (current != null) {
            path.add(0, current); // Add to the front of the list
            current = parentMap.get(current); // add to the map
        }

        // Check if the destination was reached
        if(path.size() > 0 && path.get(0) == src)
            return path;
        else    
            return new ArrayList<>();
       
    }
  
  /**
   * Finds all of the connected components in the given graph. Note
   * that for directed graphs, computes the "weakly" connected components.
   * @param g The graph to search
   * @return A mapping from vertices to components where components
   *         are numbered from 0 to n-1 for n discovered components.
   */
  public static Map<Integer,Integer> weakComponents(Graph g) {
    // @TODO: implement this function.
    // NOTE: This function CAN reuse your BFS function above.

    Map<Integer, Integer> weakComponentMap = new HashMap<>();
    boolean[] visitedVertices = new boolean[g.vertices()]; // Boolean array for visited vertices with the size equal to total number of vertices 
    int counter = 0;

        for (int i = 0; i < g.vertices(); i++) {  // look at each vertex
            if (visitedVertices[i] == false && g.hasVertex(i) == true) { // if the vertex exists and not visited yet
                // calling bfs method for vertex taht is not visited
                Map<Integer, Integer> bfsMap = bfs(g, i, false); //reusing my bfs method
                if (bfsMap != null) {
                    for (Integer vertex : bfsMap.keySet()) { //  each vertext in the bfs map
                        weakComponentMap.put(vertex, counter); //add to the weak component map
                        visitedVertices[vertex] = true; // Mark this vertex as visited
                    }
                    counter++; //increment the counter
                }
            }
        }
        return weakComponentMap; // Return the map
    }
   
  /**
   * Finds the longest shortest path in the given graph as a sequence
   * of vertices. The number of vertices returned minus one is the
   * diameter of the graph.
   * @param g The graph to search
   * @param digraph If false, treat g as undirected
   * @return The length (i.e., diameter) of the longest shortest path
   */
  public static int diameter(Graph g, boolean digraph) {
    // @TODO: implement this function
    // NOTE: You must not call your BFS function to implement diameter
    if (g.vertices() == 0) {
        return 0; // If no vertices in the graph
    }

    int maxD = 0; // initialize  diameter to 0

    // go through all vertices of the graph
    for (int startVertex = 0; startVertex < g.vertices(); startVertex++) {        
        int[] distances = new int[g.vertices()];// Array to hold distances for the current start vertex
        boolean[] visitedVertices = new boolean[g.vertices()]; // Boolean array for visited vertices

        // Initialize distances all to -1
        for (int i = 0; i < distances.length; i++) {
            distances[i] = -1; // unvisited = -1
        }
        distances[startVertex] = 0; // Distance of a vertex from itself = 0

        Queue<Integer> queue = new LinkedList<>(); // Queue for bfs
        queue.add(startVertex);
        visitedVertices[startVertex] = true;

        while (queue.isEmpty() == false) {
            int current = queue.poll();
            for (Integer neighbor : g.out(current)) { // Get outgoing neighbors
                if (visitedVertices[neighbor] == false) {
                    visitedVertices[neighbor] = true; // Make it visited
                    distances[neighbor] = distances[current] + 1; // increase distance by 1
                    queue.add(neighbor); // Add neighbor to the queue
                }
            }
            // If undirected, also check incoming edges
            if (digraph == false) {
                for (Integer neighbor : g.in(current)) { // Get incoming neighbors
                    if (visitedVertices[neighbor] == false) {
                        visitedVertices[neighbor] = true; // make it visited
                        distances[neighbor] = distances[current] + 1; // increase distance by 1
                        queue.add(neighbor); // Add neighbor to the queue
                    }
                }
            }
        }
        // Find the max of all distances
        for (int distance : distances) {
            if (distance > maxD) {
                maxD = distance;
            }
        }
    }
    return maxD; // Return the maximum distance, which is the diameter of the graph
  }  
}


