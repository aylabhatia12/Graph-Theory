/**
 * CPSC 450, Fall 2024
 * 
 * NAME: Ayla Bhatia
 * DATE: Fall 2024
 */

package cpsc450;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Basic adjacency List implementation of the Graph interface.
 */
public class AdjList implements Graph {

  private int vertexCount;                     // total number of vertices
  private int edgeCount;                       // running count of number of edges
  private Map<Integer,Set<Integer>> outEdges;  // storage for the out edges
  private Map<Integer,Set<Integer>> inEdges;   // storage for the in edges

  /**
   * Create an adjacency list (graph) given a specific (fixed) number
   * of vertices.
   * @param vertices The number of vertices of the graph.
   */
  public AdjList(int vertices) throws GraphException {
    if (vertices <= 0) {
      throw new GraphException("Number of vertices must be positive.");
    }
  
    this.vertexCount = vertices;
    this.edgeCount = 0;
    this.outEdges = new HashMap<>();
    this.inEdges = new HashMap<>();
    
    // Initialize the adjacency lists for each vertex
    for (int i = 0; i < vertices; i++) {
        outEdges.put(i, new HashSet<>());
        inEdges.put(i, new HashSet<>());
    }
  } 
  
  //Method to add an edge to list
  @Override
  public void addEdge(int x, int y) {
    // Check for valid vertex indices
    if (x >= 0 && x <  vertexCount && y >= 0 && y < vertexCount) {
      // Add the edge from x to y in the outEdges map
      if (outEdges.get(x).contains(y) == false) {
        outEdges.get(x).add(y);
        inEdges.get(y).add(x); // Add incoming edge for y
        edgeCount++;
      }     
    }  
  }
  
  //Method to remove an edge from list
  @Override
  public void removeEdge(int x, int y) {
      if (x >= 0 && x <  vertexCount && y >= 0 && y < vertexCount) {
      if (outEdges.containsKey(x) && outEdges.get(x).contains(y)) {
        outEdges.get(x).remove(y);  // Remove edge x -> y
        edgeCount--;  // Decrement edge count
    }
        
    }
  }

  //Method to get outgoing edges
  @Override
  public Set<Integer> out(int x) {
    // TODO: Implement
    Set<Integer> result = new HashSet<>();
    if (x >= 0 && x < vertexCount) {
        // Return the set of out-edges
      return new HashSet<>(outEdges.get(x));
    }
    return new HashSet<>();
  }

  //Method to get in coming edges
  @Override
  public Set<Integer> in(int x) {
    // TODO: Implement
    Set<Integer> result = new HashSet<>();
    if (x >= 0 && x < vertexCount) {
        // Return the set of in-edges
      return new HashSet<>(inEdges.get(x));
    }
    return result;
  }

  @Override
  public Set<Integer> adj(int x) {
    // TODO: Implement
    Set<Integer> result = new HashSet<>();
    if (x >= 0 && x < vertexCount) {
        Set<Integer> adjVertices = new HashSet<>(outEdges.get(x)); // Out-edges
        adjVertices.addAll(inEdges.get(x)); // In-edges
        return adjVertices;
    }
    return result;
  }

  //Method to see if an edge exists
  @Override
  public boolean hasEdge(int x, int y) {
    // TODO: Implement
    if (outEdges.containsKey(x)) {
        return outEdges.get(x).contains(y); // Return true if y is in x's adjacency list
    }
    return false;
  }

  //Method to see if a vertex exists
  @Override
  public boolean hasVertex(int x) {
    // TODO: Implement
      return (x >= 0 && x < vertexCount);
  }

  //Method to get vertex count
  @Override
  public int vertices() {
    // TODO: Implement
    return vertexCount;
  } 

  //Method to get edge count
  @Override
  public int edges() {
    // TODO: Implement
    return edgeCount;
  } 

  public void adjPrint(){
    System.out.println();
    /*for(int i=0; i<vertexCount*vertexCount; ++i){
      System.out.println()
    } */
  }
}
