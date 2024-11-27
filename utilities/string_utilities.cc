// gemc headers
#include "string_utilities.h"

// mlibrary
#include "gstring.h"

using namespace gstring;

// C++ headers
#include <algorithm>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"

using namespace CLHEP;

vector <string> get_info(string input, string chars) {
    string stripped = replaceCharInStringWithChars(input, chars, " ");

    return getStringVectorFromString(stripped);
}

vector <string> get_info(string input) {

    string chars("(),\"");

    string stripped = replaceCharInStringWithChars(input, chars, " ");

    return getStringVectorFromString(stripped);
}


// returns a vector of strings from a stringstream, space is delimeter
// ignore instances of second string
vector <string> get_strings_except(string input, string ignore) {
    vector <string> pvalues;
    stringstream plist(input);
    while (!plist.eof()) {
        string tmp;
        plist >> tmp;

        if (tmp.find(ignore) == string::npos)
            pvalues.push_back(trimSpacesFromString(tmp));
    }

    return pvalues;
}


double scan_number(const char *str) {

    if (!strcmp(str, "none")) {
        return 0;
    }

    // Scan the c_string str for numbers only, then return the value as a float.
    // The str is not allowed to have spaces or alphanumerics, only 0-9 and . but not comma
    int i = 0;

    // checking if the number is a valid digit, w/o commas (for example 1,000 is not valid). Exiting if something is wrong, except if the string is "none"
    while (char c = str[i++])
        if (!(isdigit(c) || c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E')) {
            cout << "WARNING: Unexpected alphanumeric character found in number string:" << str << endl;
            cout << "Exiting " << endl;
            exit(4);
        }

    return (stringToDouble(str));
}


// assigns G4 units to a variable
/// \fn double get_number(string v)
/// \brief Return value of the input string, which may or may not
/// contain units
/// \param v input string. Ex: 10.2*cm
/// \return value with correct G4 unit.
double get_number(string v, int warn_no_unit) {
    string value = trimSpacesFromString(v);

    if (value.find("*") == string::npos) {
        // no * found, this should be a number
        // No unit is still ok if the number is 0
        if (value.length() > 0 && warn_no_unit && stringToDouble(value) != 0) cout << "Warning: All numbers should be paired with a unit: " << v << endl;

        return scan_number(value.c_str());

    } else {
        double answer = scan_number(value.substr(0, value.find("*")).c_str());
        string units = trimSpacesFromString(value.substr(value.find("*") + 1, value.find("*") + 20));
        if (units == "m") answer *= m;
        else if (units == "inches") answer *= 2.54 * cm;
        else if (units == "inch") answer *= 2.54 * cm;
        else if (units == "cm") answer *= cm;
        else if (units == "mm") answer *= mm;
        else if (units == "um") answer *= 1E-6 * m;
        else if (units == "fm") answer *= 1E-15 * m;
        else if (units == "deg") answer *= deg;
        else if (units == "degrees") answer *= deg;
        else if (units == "arcmin") answer = answer / 60.0 * deg;
        else if (units == "rad") answer *= rad;
        else if (units == "mrad") answer *= mrad;
        else if (units == "eV") answer *= eV;
        else if (units == "MeV") answer *= MeV;
        else if (units == "KeV") answer *= 0.001 * MeV;
        else if (units == "GeV") answer *= GeV;
        else if (units == "T") answer *= tesla;
        else if (units == "T/m") answer *= tesla / m;
        else if (units == "Tesla") answer *= tesla;
        else if (units == "gauss") answer *= gauss;
        else if (units == "kilogauss") answer *= gauss * 1000;
        else if (units == "ns") answer *= ns;
        else if (units == "na") answer *= 1;
        else if (units == "counts") answer *= 1;
        else {
            cout << ">" << units << "<: unit not recognized for string <" << v << ">. Exiting" << endl;
            exit(5);

        }
        return answer;
    }

    return 0;
}


void print_vstring(vector <string> s) {
    for (unsigned int i = 0; i < s.size(); i++) {
		cout << "string element: " << i << "  content: " << s[i] << endl;
	}
}


string get_variation(string var) {
    string variation = var;
    vector <string> vars = getStringVectorFromString(var);
    if (vars[0] == "default" && vars.size() > 1) {
        variation = vars[1];
    }
    return variation;
}

bool is_main_variation(string var) {
    if (var.find("default:") == 0)
        return 1;

    return 0;
}


// overloading << for map<string, string>
ostream &operator<<(ostream &stream, map <string, string> smap) {
    cout << endl;
    for (map<string, string>::iterator it = smap.begin(); it != smap.end(); it++) {
		cout << it->first << " " << it->second << endl;
	}

    cout << endl;

    return stream;
}
