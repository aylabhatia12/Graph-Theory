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
}


