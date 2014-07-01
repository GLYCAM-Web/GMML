#ifndef TOPOLOGYBONDTYPE_HPP
#define TOPOLOGYBONDTYPE_HPP

#include <string>
#include <iostream>
#include <vector>

namespace TopologyFileSpace
{
    class TopologyBondType
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyBondType();

            TopologyBondType(int index, double force_constant, double equilibrium_value);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the bond type index
              * @return index_ attribute of the current object of this class
              */
            int GetIndex();
            /*! \fn
              * An accessor function in order to access to the bond type force constant
              * @return force_constant_ attribute of the current object of this class
              */
            double GetForceConstant();
            /*! \fn
              * An accessor function in order to access to the bond type equilibrium value
              * @return equilibrium_value_ attribute of the current object of this class
              */
            double GetEquilibriumValue();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the index of the current object
              * Set the index_ attribute of the current topology bond type
              * @param index The index attribute of the current object
              */
            void SetIndex(int index);
            /*! \fn
              * A mutator function in order to set the force constant of the current object
              * Set the force_constant_ attribute of the current topology bond type
              * @param force_constant The index attribute of the current object
              */
            void SetForceConstant(double force_constant);
            /*! \fn
              * A mutator function in order to set the force constant of the current object
              * Set the force_constant_ attribute of the current topology bond type
              * @param force_constant The force constant attribute of the current object
              */
            void SetEquilibriumValue(double equilibrium_value);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the bond type contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int index_;
            double force_constant_;
            double equilibrium_value_;

    };
}

#endif // TOPOLOGYBONDTYPE_HPP
