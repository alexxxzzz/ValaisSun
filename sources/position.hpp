#ifndef POSITION_HEADER_DEFINED
#define POSITION_HEADER_DEFINED

struct position {
    uint32_t x;
    uint32_t y;
    position() {}
    position(uint32_t x_, uint32_t y_) : x(x_), y(y_) {}
    bool operator<(const position& p) const
    {
        uint64_t i1;
        uint64_t i2;
        i1 = x + ((uint64_t)y << 32);
        i2 = p.x + ((uint64_t)p.y << 32);
        return i1 < i2;
    }
    bool operator==(const position& p) const
    {
        if (x != p.x) return false;
        if (y != p.y) return false;
        return true;
    }
};

#endif //POSITION_HEADER_DEFINED