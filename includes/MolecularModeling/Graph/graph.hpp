#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "node.hpp"
#include "edge.hpp"

namespace GraphDS
{
 template < class N > class Graph
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
             * List of nodes in a Graph
             */
            typedef typename std::vector<GraphDS::Node<N>*> NodeVector;
            typedef typename std::vector<GraphDS::Edge<N>*> EdgeVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
             Graph();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////

             /*! \fn
               * An accessor function in order to access the nodes of the graph
               * @return  graphNodeList_ attribute of the current object of Graph class
               */
             NodeVector GetGraphNodeList();

             /*! \fn
               * An accessor function in order to access the edges of the graph
               * @return  graphEdgeList_ attribute of the current object of Graph class
               */
             EdgeVector GetGraphEdgeList();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////

             /*! \fn
               * A mutator function in order to set nodes to the current Graph
               * @param graphNodeList_ List of nodes of the current Graph class
               */
             void SetGraphNodeList(NodeVector nodeList);

             /*! \fn
               * A mutator function in order to set edges to the current Graph
               * @param graphEdgeList_ List of outgoing/incoming edges of the current Graph class
               */
             void SetGraphEdgeList(EdgeVector edgeList);

            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to add a new node to the Graph
              * Set the node_ attributes of the current graph object
              * @param node_ The node attribute of the current graph object
              */
           void AddNewNode(GraphDS::Node<N>* newNode);

           /*! \fn
             * A function to remove an existing new node from the Graph
             * Remove the node_ attributes of the current graph object
             * @param node_ The node attribute of the current graph object
             */
          void RemoveNode(GraphDS::Node<N>* delNode);

           /*! \fn
             * A function to add an edge between two nodes of the current graph object of this class
             * @param node_ The node attribute of the current graph object
             */
          void AddEdge(GraphDS::Node<N>* firstNode, GraphDS::Node<N>* secondNode);

           /*! \fn
             * A function to find a node in the current graph object using the node id
             * @return node_ attribute of the current graph object of this class
             */

           GraphDS::Node<N>* FindNodeById(std::string node_id);

            /////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the nodes in a graph.
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            NodeVector graphNodeList_;         /*!< List of all nodes in the graph >*/
            EdgeVector graphEdgeList_;         /*!< List of all edges in the graph >*/
    };
}
#include "../../../src/MolecularModeling/Graph/graph.cc"
#endif // GRAPH_HPP
