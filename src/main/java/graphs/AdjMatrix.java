/**
 * 
 * NAME: Ayla Bhatia
 */

package graphs;

import java.util.HashSet;
import java.util.Set;

/**
 * Adjacency Matrix implementation of the Graph interface. 
 */
public class AdjMatrix implements Graph {

  private int vertexCount;      // total number of vertices
  private int edgeCount;        // running count of number of edges
  private boolean matrix[];     // storage for the matrix as "flattened" 2D array

  /**
   * Create an adjacency matrix (graph) given a specific (fixed)
   * number of vertices.
   * @param vertices The number of vertices in the graph. 
   */ 
  public AdjMatrix(int vertices) throws GraphException {
    if (vertices <= 0) {
      throw new GraphException("Number of vertices must be positive.");
    }
  
    this.vertexCount = vertices;
    this.edgeCount = 0;
    this.matrix = new boolean[vertices * vertices]; // Initialize flattened 2D array  
  }

  //New method//
  public void adjPrint(){
    
    for(int i=0; i<vertexCount*vertexCount; ++i){
      System.out.print(matrix[i]+ " " );
    } 
     System.out.println();
  }

  //Method to add an edge to the matrix
  @Override
  public void addEdge(int x, int y) {
    if (x >= 0 && x < vertexCount && y >= 0 && y < vertexCount) {
      if (!matrix[vertexCount * x + y] ) {    
        matrix[vertexCount*x+y]=true;
       // matrix[vertexCount*y+x]=true;
        edgeCount++;
      }
     // if (!matrix[vertexCount * y + x]) {    
     //   matrix[vertexCount*y+x]=true;
      //  edgeCount++;
     // }
    }
  }

  //Method to remove an edge from the matrix
  @Override
  public void removeEdge(int x, int y) {
   if (x >= 0 && x < vertexCount && y >= 0 && y < vertexCount) {     
      if (matrix[x * vertexCount + y]) {
          matrix[x * vertexCount + y] = false;
          edgeCount--;
      }
      //if (matrix[y * vertexCount + x]) {
      //    matrix[y * vertexCount + x] = false;
      //    edgeCount--;
      //}
    }
  }


  // Method to get outgoing edges from the matrix
  @Override
  public Set<Integer> out(int x) {
    Set<Integer> result = new HashSet<>();
    for (int i = 0; i < vertexCount; i++) {
        if (matrix[x * vertexCount + i]) {
            result.add(i);
        }
    }
    return result;
  }

  // Method to get incoming edges from the matrix
  @Override
  public Set<Integer> in(int x) {
    Set<Integer> result = new HashSet<>();
    for (int i = 0; i < vertexCount; i++) {
        if (matrix[i * vertexCount + x]) {
            result.add(i);
        }
    }
    return result;
  }

  // Method to get adjacent edges from the matrix
  @Override
 public Set<Integer> adj(int x) {
    Set<Integer> result = new HashSet<>();
    // Check outgoing edges
    for (int i = 0; i < vertexCount; i++) {
        if (matrix[x * vertexCount + i] == true ||  matrix[i * vertexCount + x] == true) {
            result.add(i);
        }
    }
    // Check incoming edges
   /* for (int i = 0; i < vertexCount; i++) {
        if (matrix[i * vertexCount + x]) {
            result.add(i);
        }
    }*/
    return result;
  }

  //Method to check if an edge exists
  @Override
  public boolean hasEdge(int x, int y) {
    if (x >= 0 && x < vertexCount && y >= 0 && y < vertexCount){
      return matrix[vertexCount*x+y];
    }
    else{
      return false;
    }
  }

  //Method to see if a vertex exists
  @Override
  public boolean hasVertex(int x) {
    return (x >= 0 && x < vertexCount);
  }

  //Method to get vertext count
  @Override
  public int vertices() {
    return vertexCount;
  }
  
  //Method to get edge count
  @Override
  public int edges() {
    return edgeCount;
  }  
}
