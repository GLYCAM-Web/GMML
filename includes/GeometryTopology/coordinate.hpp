#ifndef COORDINATE_HPP
#define COORDINATE_HPP

#include <ostream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>

namespace GeometryTopology
{
    class Coordinate;
    typedef std::vector<Coordinate*> CoordinateVector;
    class Coordinate
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Coordinate();
            /*! \fn
              * Constructor with given three double values
              * @param x A double value on X axis in cartesian coordinate
              * @param y A double value on Y axis in cartesian coordinate
              * @param z A double value on Z axis in cartesian coordinate
              */
            Coordinate(const std::string x, const std::string y, const std::string z);
            Coordinate(double x, double y, double z);
            /*! \fn
              * Copy constructor
              * @param coordinate A coordinate to be copied to another instance
              */
            Coordinate(const Coordinate& coordinate);
            Coordinate(Coordinate* coordinate);


            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to X coordinate attribute of the current object
              * The attribute is set by the contents of the given file
              * @return x_ of the current object of this class
              */
            double GetX() const;
            /*! \fn
              * An accessor function in order to access to Y coordinate attribute of the current object
              * The attribute is set by the contents of the given file
              * @return y_ of the current object of this class
              */
            double GetY() const;
            /*! \fn
              * An accessor function in order to access to Z coordinate attribute of the current object
              * The attribute is set by the contents of the given file
              * @return z_ of the current object of this class
              */
            double GetZ() const;
/** @}*/
            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            void SetX(double x);
            void SetY(double y);
            void SetZ(double z);
/** @}*/
            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Verifiers_and_Issue_Resolvers
               * @{
               */
            /*! \fn
              * Translate a coordinate with respect to given coordinate by x, y and z
              * @param x A double value on X axis for the origin of translation
              * @param y A double value on Y axis for the origin of translation
              * @param z A double value on Z axis for the origin of translation
              */
            void Translate(double x, double y, double z);
            /*! \fn
              * Compare current object with the given coordinate and return true if they are the same
              * @param coordinate A coordinate to be compared with the current object
              * @return A boolean value as a result of comparison
              */
            bool CompareTo(Coordinate coordinate) const;
            bool withinDistance(const Coordinate *coordinate, const double distance) const;

            /*! \fn
              * Calculate the distance between current object coordinate and the given coordinate
              * @param coordinate A coordinate to calculate the distance between that and the current object
              * @return A double value as a distance between two coordinates
              */
            double Distance(const Coordinate &coordinate) const;

            /*! \fn
              * Calculate the distance between current object coordinate and the given coordinate
              * @param coordinate A coordinate to calculate the distance between that and the current object
              * @return A double value as a distance between two coordinates
              */
            double Distance(const Coordinate *coordinate) const;

            /*! \fn
              * Calculate the length of the vector
              * @return Length of the vector
              */
            double length() const;
            /*! \fn
              * Normalize the vector and update it on the current object
              */
            void Normalize();
            /*! \fn
              * Calculate the dot product of current object and the given one
              * @param coordinate Second vector in the dot product operator
              * @return Dot product value of two vectors
              */
            double DotProduct(Coordinate coordinate);
            /*! \fn
              * Calculate the vector from cross producting the current object and the given one
              * @param coordinate Second vector in the cross product operator
              * @return Cross product vector of two vectors
              */
            void CrossProduct(Coordinate coordinate);
            /*! \fn
              * Add current object to the given one and update the current one
              * @param coordinate Second vector in the add operator
              */
/** @}*/
            void operator+(Coordinate coordinate);
            void operator +(double addition);
            /*! \fn
              * Subtract current object to the given one and update the current one
              * @param coordinate Second vector in the subtract operator
              */
            void operator-(Coordinate coordinate);
            void operator /(Coordinate coordinate);
            void operator/(double divisor);
            void operator*(double multiplier);
            void TranslateAll(CoordinateVector coordinate_set, double margin, int pos);
            void RotateAngularAll(CoordinateVector coordinate_set, double angle, int pos);
            void RotateTorsionalAll(CoordinateVector coordinate_set, double torsion, int pos);

            /*! \fn
             * Convert internal coordinate to the corresponding cartesian coordinate
             * @param coordinate_list List of at most three internal coordinates in order to calculate 
             * the cartesian coordinate of the given internal coordinate (distance, angle, torsion)
             * @param distance X value of the internal coordinate
             * @param angle Y value of the interanl coordinate
             * @param torsion Z value of the internal coordinate
             * @return Cartesian coordinate of the internal coordinate (distance, angle, torsion)
             */
            Coordinate* ConvertInternalCoordinate2CartesianCoordinate(
                       CoordinateVector coordinate_list, 
                       double distance, double angle, double torsion);
            Coordinate* ConvertCartesianCoordinate2InternalCoordinate(
                       Coordinate* coordinate, CoordinateVector coordinate_list);
            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cerr) const;
            std::string ToString() const;
	private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            double x_;          /*!< x */
            double y_;          /*!< y */
            double z_;          /*!< z */
    };
}

#endif // COORDINATE_HPP
