#ifdef NDEBUG
#define DEBUG(x)
#else
#define DEBUG(x) do {std::cerr << x << std::endl << std::flush; } while(0)
#endif