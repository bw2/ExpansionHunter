#include <string>


std::string add_commas_at_thousands(std::string n) {
    if (n.length() > 3) {
        // Loop to insert commas for thousands separators in the string
        for (int i = n.length() - 3; i > 0; i -= 3) {
            n.insert(i, ",");
        }
    }

    return n;
}

std::string add_commas_at_thousands(int n) {
	return add_commas_at_thousands(std::to_string(n));
}

std::string add_commas_at_thousands(long n) {
	return add_commas_at_thousands(std::to_string(n));
}

std::string add_commas_at_thousands(unsigned long n) {
	return add_commas_at_thousands(std::to_string(n));
}
