#ifndef p2pt_h_
#define p2pt_h_

namespace P2Pt {

// At most 8 solutions with positive depth, TODO: assert if longer
static constexpr int TS_MAX_LEN = 8;
static constexpr int RT_MAX_LEN = (TS_MAX_LEN * TS_MAX_LEN);

template <typename T=double>
class p2pt { // fully static, not to be instantiated - just used for templating
	public:
  static void pose_from_point_tangents(
    const T (&gama1)[3], const T (&tgt1)[3],
    const T (&gama2)[3], const T (&tgt2)[3],
    const T (&Gama1)[3], const T (&Tgt1)[3],
    const T (&Gama2)[3], const T (&Tgt2)[3],
    T (*output_RT)[RT_MAX_LEN][4][3],
    int *output_RT_len,
    T *output_degen
  );
};

} // namespace P2Pt

#endif  // !p2pt_h_

