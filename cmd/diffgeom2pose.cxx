#include "p2pt/p2pt.h"

using namespace P2Pt;

#define Float double
static constexpr Float tol = 1e-3;


// Exit code. Conventions:

void
print_usage()
{
  std::cerr << "Usage: diffgeom2pose input solutions\n\n";
  std::cerr << "If no argument is given, 'input' is assumed stdin,\n\
  solutions will be output to stdout\n";
  std::cerr << "Example: \n"
               "  diffgeom2pose input_file solutions_file\n"
               "  minus <input_file >solutions_file\n"
               "  minus -g       # (or --profile) : performs a default solve for profiling\n"
               "  minus -i       # (or --image_data) : reads point-tangents from stdin\n"
               "  minus -h       # (or --help) : print this help message\n"
            <<
  R"(-i | --image_data usage:
 
  Input format (notation _points_coords. any number of spaces and newlines optional. can be in
  one row or one column as well).  
 
  p00 p01  
  p10 p11 
  
  t00 t01
  t10 t11

  P00 P01 P02 
  P10 P11 P12

  T00 T01 T02 
  T10 T11 T12
  
  K00 K01 K02       # intrinsic parameters: only these elements
   0  K11 K22

  r000 r001 r002    # GROUND TRUTH (optional) if -gt flag provided, pass the ground truth here:
  r010 r011 r012    # default camera format if synthcurves flag passed: 
  r020 r021 r022    # just like a 3x4 [R|T] but transposed to better fit row-major:
   c00  c01  c02    #         | R |
                    # P_4x3 = | - |
                    #         | C'|

  # One way to use this is 
  #     synthdata | minus-chicago -i
  # where synthdata is provided in diffgeom2pose/scripts)";
             
  exit(1);
}

bool stdio_ = true;  // by default read/write from stdio
bool ground_truth_ = false;
std::ifstream infp_;
bool image_data_ = false;
const char *input_ = "stdin";
const char *output_ = "stdout";


// Output solutions in ASCII matlab format
//
// ---------------------------------------------------------
// If in the future our solver is really fast, we may need Binary IO:
// complex solutions[NSOLS*NVE];
// 
// To read this output file in matlab, do:
// fid = fopen(fname,'r');
// a_raw = fread(fid,'double');
// fclose(fid);
//
// Reshape a to have proper real and imaginary parts
// a = a_raw(1:2:end) + i*a_raw(2:2:end);
// 
template <typename F=double>
static bool
mwrite(const M::solution s[M::nsols], const char *fname)
{
  bool scilab=false;
  std::string imag("+i*");
  if (scilab) imag = std::string("+%i*");
    
  std::ofstream fsols;
  std::streambuf *buf;
  
  if (stdio_) {
    buf = std::cout.rdbuf();
    std::cout << std::setprecision(20);
  } else {
    fsols.open(fname,std::ios::out);
    if (!fsols) {
      std::cerr << "minus: error, unable to open file name" << std::endl;
      return false;
    }
    buf = fsols.rdbuf();
    fsols << std::setprecision(20);
  }
  
  std::ostream out(buf);
  out << std::setprecision(20);
  out << "[";
  for (unsigned i=0; i <M::nsols; ++i) {
    for (unsigned var=0; var < M::nve; ++var) {
      out << s[i].x[var].real() << imag << s[i].x[var].imag();
      if (i*var +1 < M::nve * M::nsols) 
        out << std::endl;
      // BINARY fsols.write((char *)(s[i].x[var]),2*sizeof(double));
    }
  }
  out << "]\n";
  
  if (!stdio_) fsols.close();
  return true;
}
// Try to read n elements, filling in p in row-major order.
template <typename F=double>
static bool
read_block(std::istream &in, F *p, unsigned n)
{
  LOG("reading");
  const F *end = p + n;
  while (!in.eof() && p != end) {
      try {
        in >> *p++;
//        std::cerr << *(p-1) << std::endl;
        if (in.eof()) {
          std::cerr << "I/O Error: Premature input termination\n";
          return false;
        }
      } catch (std::istream::failure &E) {
        std::cerr << "I/O Error: Invalid input conversion or other error\n";
        return false;
      }
  }
  if (p != end) {
    std::cerr << "I/O Premature input termination\n";
    return false;
  }
  return true;
}

static bool
init_input(const char *fname, std::istream *inp)
{
  if (!stdio_) {
    infp_.open(fname, std::ios::in);
    if (!infp_) {
      std::cerr << "I/O Error opening input " << fname << std::endl;
      return false;
    }
    inp = &infp_;
  }
  inp->exceptions(std::istream::failbit | std::istream::badbit);
  return true;
}
  
//
// Reads the format specified in the print_usage() for the -i flag
// 
// This is processed into the global params_start_target_
// 
template <typename F=double>
static bool
iread(std::istream &in)
{
  LOG("reading p_");
  if (!read_block(in, (F *)data::p_, io::pp::nviews*io::pp::npoints*io::ncoords2d))
    return false;
  LOG("reading tgt_");
  if (!read_block(in, (F *)data::tgt_, io::pp::nviews*io::pp::npoints*io::ncoords2d))
    return false;
  unsigned tgt_ids[2];
  LOG("reading tgt_ids");
  if (!read_block<unsigned>(in, tgt_ids, 2))
    return false;
  if (reading_first_point_) {
    LOG("reading K_");
    if (!read_block(in, (F *) data::K_, io::ncoords2d*io::ncoords2d_h))
      return false;
    LOG("reading ground truth cams");
    if (ground_truth_ && !read_block(in, (F *) data::cameras_gt_, io::pp::nviews*4*3))
      return false;
    io::point_tangents2params_img(data::p_, data::tgt_, tgt_ids[0], tgt_ids[1],
        data::K_, data::params_start_target_);
    reading_first_point_ = false;
  } else { // when reading second point B, do not gammify A again
    static constexpr bool gammify_target_problem = false;
    io::point_tangents2params_img(data::p_, data::tgt_, tgt_ids[0], tgt_ids[1],
        data::K_, data::params_start_target_, gammify_target_problem);
  }
  return true;
}

// reads into the global variable params_
// Format is just like P01 variable in solveChicago in chicago.m2
// and contains the concatenated parameters of the start system
// and of the target system, with some randomization to improve conditioning.
// But here there is no imaginary 'i' string:
//
// P01(0).real()  P01(0).imag()      // I mean P01(0) or P01#0
// P01(1).real()  P01(1).imag()
// P01(2).real()  P01(2).imag()
// P01(3).real()  P01(3).imag()
// ...
//
// The file can also be one line, listing the above in row-major order like so:
// 
// P01(0).real()  
// P01(0).imag() 
// P01(1).real()
// P01(1).imag()
// ...
// 
// It is up to the user to build this from an actual input for a target system,
// be it point-tangents as in Ric's format, be it a linecomplex as in Hongy's format
//
// This format is generic enough to be adapted to M2 or matlab
template <typename F=double>
static bool
mread(std::istream &in)
{
  F *dparams = (F *)data::params_;
  while (!in.eof() && dparams != (F *)data::params_+2*2*M::f::nparams) {
      try {
      in >> *dparams++;
      // std::cerr << "reading " <<  *(dparams-1) << std::endl;;
      if (in.eof()) {
        std::cerr << "I/O Error: Premature input termination\n";
        return false;
      }
      in >> *dparams++;
      } catch (std::istream::failure &E) {
        std::cerr << "I/O Error: Invalid input conversion or other error\n";
        return false;
      }
  }
  if (dparams != (F *)data::params_+2*2*M::f::nparams)
    std::cerr << "I/O Premature input termination\n";
//  for (unsigned i=0; i < 2*NPARAMS; ++i)
//    std::cerr << "D " << params_[i] << std::endl;
  return true;
}

void
process_args(int argc, char **argv)
{
  --argc; ++argv;
  // switches that can show up only in 1st position
  
  enum {INITIAL_ARGS, AFTER_INITIAL_ARGS, IMAGE_DATA} argstate = INITIAL_ARGS;
  bool incomplete = false;
  std::string arg;
  if (argc) {
    arg = std::string(*argv);
    if (arg == "-h" || arg == "--help")
      print_usage();
    } else if (arg == "-i" || arg == "--image_data") {
      image_data_ = true; 
      argstate = IMAGE_DATA;
      --argc; ++argv;
    } else if (arg[0] != '-') {
      if (argc == 2) {
          input_ = argv[1];
          output_ = argv[2];
          stdio_ = false;
      } else {
          std::cerr << "minus: \033[1;91m error\e[m\n";
          print_usage();
      }
    }
    
    while (argc) { // second and beyond: above switches must already be set
      arg = std::string(*argv);
      LOG("parsing arg " + arg);
      
      // argstate >= AFTER_INITIAL_ARGS ----------------------------------------
      if (argstate == IMAGE_DATA) {
        if (arg == "-gt") {
          ground_truth_ = true;
          --argc; ++argv;
          argstate = IMAGE_DATA;
          continue;
        }
        argstate = AFTER_INITIAL_ARGS;
        continue;
      }
      
      // argstate == AFTER_INITIAL_ARGS ----------------------------------------
      std::cerr << "minus: \033[1;91m error\e[m\n - unrecognized argument " << arg << std::endl;;
      print_usage();
    }

  }
}

// Simplest possible command to compute the Chicago problem
// for estimating calibrated trifocal geometry from points and lines at points
//
// This is to be kept very simple C with only minimal C++ with Templates.
// If you want to complicate this, please create another executable.
// 
int
main(int argc, char **argv)
{
  std::istream *inp = &std::cin;
  
  process_args(argc, argv);

  if (image_data_) {
    LOG("param: input is image pixel data");
    if (ground_truth_)
      LOG("param: reading ground truth appended to input pixel data");
  }
  if (stdio_)
    LOG("reading from stdio");
  else
    LOG("reading from " << input_ << " writing to " << output_);

  init_input(input_, inp);
  if (image_data_) {  // read image pixel-based I/O parameters
    if (!iread<Float>(*inp))
      return 1;
  } else {  // read raw I/O homotopy parameters (to be used as engine)
    if (!mread<Float>(*inp))  // reads into global params_
      return 1;
  }
  
  {
    std::cerr << "LOG \033[0;33mStarting solver\e[m\n" << std::endl;

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    // XXX run solver
     
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t2 - t1).count();
    std::cerr << "LOG \033[1;32mTime of solver: " << duration << "ms\e[m" << std::endl;
  }

  if (!mwrite<Float>(solutions, output_)) return 2;

  // ---------------------------------------------------------------------------
  // test_final_solve_against_ground_truth(solutions);
  // optional: filter solutions using positive depth, etc.
  if (ground_truth_) {
    // TODO(juliana) should we has_valid_solutions here? 
    io::RC_to_QT_format(data::cameras_gt_, data::cameras_gt_quat_);
    unsigned sol_id;
    bool found = io::probe_all_solutions(solutions, data::cameras_gt_quat_, &sol_id);
    if (found) {
      LOG("found solution at index: " << sol_id);
      if (solutions[sol_id].status != M::REGULAR)
        LOG("PROBLEM found ground truth but it is not REGULAR: " << sol_id);
    } else {
      LOG("\033[1;91mFAIL:\e[m  ground-truth not found among solutions");
      return SOLVER_FAILURE; 
      // you can detect solver failure by checking this exit code.
      // if you use shell, see:
      // https://www.thegeekstuff.com/2010/03/bash-shell-exit-status
    }
  } else if (!has_valid_solutions(solutions)) { // if no ground-truth is provided, it will return error
    LOG("\033[1;91mFAIL:\e[m  no valid solutions");
    return SOLVER_FAILURE;                    // if it can detect that the solver failed by generic tests
  }
  return 0;
}
