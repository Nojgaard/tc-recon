#ifndef TC_RANGE_HPP
#define TC_RANGE_HPP
#include <utility>

namespace tc {

template <typename Iter>
class Range {
public:
	Range(std::pair<Iter, Iter> pi): _begin(pi.first), _end(pi.second) {}

	Iter begin() const { return _begin; }
	Iter end() const { return _end; }

private:
	Iter _begin, _end;
};

template <typename Iter>
Range<Iter> make_range(std::pair<Iter,Iter> pi) {
	return Range<Iter>(pi);
}

} /* namespace tc */
#endif /* ifndef TC_RANGE_HPP */
