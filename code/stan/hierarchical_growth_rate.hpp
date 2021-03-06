
// Code generated by stanc v2.26.1
#include <stan/model/model_header.hpp>
namespace hierarchical_growth_rate_model_namespace {


inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}

inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}


using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math;
using stan::math::pow; 

stan::math::profile_map profiles__;
static int current_statement__= 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 14, column 4 to column 31)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 15, column 4 to column 27)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 18, column 4 to column 29)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 19, column 4 to column 32)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 20, column 4 to column 25)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 21, column 4 to column 28)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 24, column 4 to column 24)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 29, column 4 to column 28)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 30, column 4 to column 26)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 31, column 4 to column 32)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 32, column 4 to column 32)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 33, column 4 to column 36)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 36, column 4 to column 36)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 37, column 4 to column 48)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 40, column 4 to column 62)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 3, column 4 to column 19)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 4, column 4 to column 19)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 5, column 30 to column 31)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 5, column 4 to column 33)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 8, column 20 to column 21)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 8, column 4 to column 26)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 9, column 20 to column 21)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 9, column 4 to column 28)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 14, column 20 to column 21)",
                                                      " (in '/Users/gchure/Dropbox/git/postdoc_projects/ltee_diauxie/code/stan/hierarchical_growth_rate.stan', line 15, column 20 to column 21)"};



class hierarchical_growth_rate_model final : public model_base_crtp<hierarchical_growth_rate_model> {

 private:
  int J;
  int N;
  std::vector<int> idx;
  Eigen::Matrix<double, -1, 1> OD;
  Eigen::Matrix<double, -1, 1> time;
 
 public:
  ~hierarchical_growth_rate_model() { }
  
  inline std::string model_name() const final { return "hierarchical_growth_rate_model"; }

  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.26.1", "stancflags = "};
  }
  
  
  hierarchical_growth_rate_model(stan::io::var_context& context__,
                                 unsigned int random_seed__ = 0,
                                 std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "hierarchical_growth_rate_model_namespace::hierarchical_growth_rate_model";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 16;
      context__.validate_dims("data initialization","J","int",
          context__.to_vec());
      J = std::numeric_limits<int>::min();
      
      current_statement__ = 16;
      J = context__.vals_i("J")[(1 - 1)];
      current_statement__ = 16;
      current_statement__ = 16;
      check_greater_or_equal(function__, "J", J, 1);
      current_statement__ = 17;
      context__.validate_dims("data initialization","N","int",
          context__.to_vec());
      N = std::numeric_limits<int>::min();
      
      current_statement__ = 17;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 17;
      current_statement__ = 17;
      check_greater_or_equal(function__, "N", N, 1);
      current_statement__ = 18;
      validate_non_negative_index("idx", "N", N);
      current_statement__ = 19;
      context__.validate_dims("data initialization","idx","int",
          context__.to_vec(N));
      idx = std::vector<int>(N, std::numeric_limits<int>::min());
      
      current_statement__ = 19;
      assign(idx, nil_index_list(), context__.vals_i("idx"),
        "assigning variable idx");
      current_statement__ = 19;
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        current_statement__ = 19;
        current_statement__ = 19;
        check_greater_or_equal(function__, "idx[sym1__]", idx[(sym1__ - 1)],
                               1);}
      current_statement__ = 19;
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        current_statement__ = 19;
        current_statement__ = 19;
        check_less_or_equal(function__, "idx[sym1__]", idx[(sym1__ - 1)], J);
      }
      current_statement__ = 20;
      validate_non_negative_index("OD", "N", N);
      current_statement__ = 21;
      context__.validate_dims("data initialization","OD","double",
          context__.to_vec(N));
      OD = Eigen::Matrix<double, -1, 1>(N);
      stan::math::fill(OD, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> OD_flat__;
        current_statement__ = 21;
        assign(OD_flat__, nil_index_list(), context__.vals_r("OD"),
          "assigning variable OD_flat__");
        current_statement__ = 21;
        pos__ = 1;
        current_statement__ = 21;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 21;
          assign(OD, cons_list(index_uni(sym1__), nil_index_list()),
            OD_flat__[(pos__ - 1)], "assigning variable OD");
          current_statement__ = 21;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 21;
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        current_statement__ = 21;
        current_statement__ = 21;
        check_greater_or_equal(function__, "OD[sym1__]", OD[(sym1__ - 1)], 0);
      }
      current_statement__ = 22;
      validate_non_negative_index("time", "N", N);
      current_statement__ = 23;
      context__.validate_dims("data initialization","time","double",
          context__.to_vec(N));
      time = Eigen::Matrix<double, -1, 1>(N);
      stan::math::fill(time, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> time_flat__;
        current_statement__ = 23;
        assign(time_flat__, nil_index_list(), context__.vals_r("time"),
          "assigning variable time_flat__");
        current_statement__ = 23;
        pos__ = 1;
        current_statement__ = 23;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 23;
          assign(time, cons_list(index_uni(sym1__), nil_index_list()),
            time_flat__[(pos__ - 1)], "assigning variable time");
          current_statement__ = 23;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 23;
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        current_statement__ = 23;
        current_statement__ = 23;
        check_greater_or_equal(function__, "time[sym1__]",
                               time[(sym1__ - 1)], 0);}
      current_statement__ = 24;
      validate_non_negative_index("OD_init", "J", J);
      current_statement__ = 25;
      validate_non_negative_index("lam", "J", J);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += J;
      num_params_r__ += J;
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += 1;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI, stan::require_vector_like_t<VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "hierarchical_growth_rate_model_namespace::log_prob";
(void) function__;  // suppress unused var warning

    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning

    
    try {
      Eigen::Matrix<local_scalar_t__, -1, 1> OD_init;
      OD_init = Eigen::Matrix<local_scalar_t__, -1, 1>(J);
      stan::math::fill(OD_init, DUMMY_VAR__);
      
      current_statement__ = 1;
      OD_init = in__.vector(J);
      current_statement__ = 1;
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        current_statement__ = 1;
        if (jacobian__) {
          current_statement__ = 1;
          assign(OD_init, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(OD_init[(sym1__ - 1)], 0, lp__),
            "assigning variable OD_init");
        } else {
          current_statement__ = 1;
          assign(OD_init, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(OD_init[(sym1__ - 1)], 0),
            "assigning variable OD_init");
        }}
      Eigen::Matrix<local_scalar_t__, -1, 1> lam;
      lam = Eigen::Matrix<local_scalar_t__, -1, 1>(J);
      stan::math::fill(lam, DUMMY_VAR__);
      
      current_statement__ = 2;
      lam = in__.vector(J);
      current_statement__ = 2;
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        current_statement__ = 2;
        if (jacobian__) {
          current_statement__ = 2;
          assign(lam, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(lam[(sym1__ - 1)], 0, lp__),
            "assigning variable lam");
        } else {
          current_statement__ = 2;
          assign(lam, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(lam[(sym1__ - 1)], 0),
            "assigning variable lam");
        }}
      local_scalar_t__ OD_init_mu;
      OD_init_mu = DUMMY_VAR__;
      
      current_statement__ = 3;
      OD_init_mu = in__.scalar();
      current_statement__ = 3;
      if (jacobian__) {
        current_statement__ = 3;
        OD_init_mu = stan::math::lb_constrain(OD_init_mu, 0, lp__);
      } else {
        current_statement__ = 3;
        OD_init_mu = stan::math::lb_constrain(OD_init_mu, 0);
      }
      local_scalar_t__ OD_init_sigma;
      OD_init_sigma = DUMMY_VAR__;
      
      current_statement__ = 4;
      OD_init_sigma = in__.scalar();
      current_statement__ = 4;
      if (jacobian__) {
        current_statement__ = 4;
        OD_init_sigma = stan::math::lb_constrain(OD_init_sigma, 0, lp__);
      } else {
        current_statement__ = 4;
        OD_init_sigma = stan::math::lb_constrain(OD_init_sigma, 0);
      }
      local_scalar_t__ lam_mu;
      lam_mu = DUMMY_VAR__;
      
      current_statement__ = 5;
      lam_mu = in__.scalar();
      current_statement__ = 5;
      if (jacobian__) {
        current_statement__ = 5;
        lam_mu = stan::math::lb_constrain(lam_mu, 0, lp__);
      } else {
        current_statement__ = 5;
        lam_mu = stan::math::lb_constrain(lam_mu, 0);
      }
      local_scalar_t__ lam_sigma;
      lam_sigma = DUMMY_VAR__;
      
      current_statement__ = 6;
      lam_sigma = in__.scalar();
      current_statement__ = 6;
      if (jacobian__) {
        current_statement__ = 6;
        lam_sigma = stan::math::lb_constrain(lam_sigma, 0, lp__);
      } else {
        current_statement__ = 6;
        lam_sigma = stan::math::lb_constrain(lam_sigma, 0);
      }
      local_scalar_t__ sigma;
      sigma = DUMMY_VAR__;
      
      current_statement__ = 7;
      sigma = in__.scalar();
      current_statement__ = 7;
      if (jacobian__) {
        current_statement__ = 7;
        sigma = stan::math::lb_constrain(sigma, 0, lp__);
      } else {
        current_statement__ = 7;
        sigma = stan::math::lb_constrain(sigma, 0);
      }
      {
        current_statement__ = 8;
        lp_accum__.add(normal_lpdf<propto__>(sigma, 0, 0.01));
        current_statement__ = 9;
        lp_accum__.add(std_normal_lpdf<propto__>(lam_mu));
        current_statement__ = 10;
        lp_accum__.add(normal_lpdf<propto__>(lam_sigma, 0, 0.01));
        current_statement__ = 11;
        lp_accum__.add(normal_lpdf<propto__>(OD_init_mu, 0, 0.1));
        current_statement__ = 12;
        lp_accum__.add(normal_lpdf<propto__>(OD_init_sigma, 0, 0.01));
        current_statement__ = 13;
        lp_accum__.add(normal_lpdf<propto__>(lam, lam_mu, lam_sigma));
        current_statement__ = 14;
        lp_accum__.add(
          normal_lpdf<propto__>(OD_init, OD_init_mu, OD_init_sigma));
        current_statement__ = 15;
        lp_accum__.add(
          normal_lpdf<propto__>(OD,
            elt_multiply(
              rvalue(OD_init, cons_list(index_multi(idx), nil_index_list()),
                "OD_init"),
              stan::math::exp(
                elt_multiply(time,
                  rvalue(lam, cons_list(index_multi(idx), nil_index_list()),
                    "lam")))), sigma));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr>
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "hierarchical_growth_rate_model_namespace::write_array";
(void) function__;  // suppress unused var warning

    (void) function__;  // suppress unused var warning

    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning

    
    try {
      Eigen::Matrix<double, -1, 1> OD_init;
      OD_init = Eigen::Matrix<double, -1, 1>(J);
      stan::math::fill(OD_init, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 1;
      OD_init = in__.vector(J);
      current_statement__ = 1;
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        current_statement__ = 1;
        assign(OD_init, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_constrain(OD_init[(sym1__ - 1)], 0),
          "assigning variable OD_init");}
      Eigen::Matrix<double, -1, 1> lam;
      lam = Eigen::Matrix<double, -1, 1>(J);
      stan::math::fill(lam, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 2;
      lam = in__.vector(J);
      current_statement__ = 2;
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        current_statement__ = 2;
        assign(lam, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_constrain(lam[(sym1__ - 1)], 0),
          "assigning variable lam");}
      double OD_init_mu;
      OD_init_mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      OD_init_mu = in__.scalar();
      current_statement__ = 3;
      OD_init_mu = stan::math::lb_constrain(OD_init_mu, 0);
      double OD_init_sigma;
      OD_init_sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      OD_init_sigma = in__.scalar();
      current_statement__ = 4;
      OD_init_sigma = stan::math::lb_constrain(OD_init_sigma, 0);
      double lam_mu;
      lam_mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 5;
      lam_mu = in__.scalar();
      current_statement__ = 5;
      lam_mu = stan::math::lb_constrain(lam_mu, 0);
      double lam_sigma;
      lam_sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 6;
      lam_sigma = in__.scalar();
      current_statement__ = 6;
      lam_sigma = stan::math::lb_constrain(lam_sigma, 0);
      double sigma;
      sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 7;
      sigma = in__.scalar();
      current_statement__ = 7;
      sigma = stan::math::lb_constrain(sigma, 0);
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        vars__.emplace_back(OD_init[(sym1__ - 1)]);}
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        vars__.emplace_back(lam[(sym1__ - 1)]);}
      vars__.emplace_back(OD_init_mu);
      vars__.emplace_back(OD_init_sigma);
      vars__.emplace_back(lam_mu);
      vars__.emplace_back(lam_sigma);
      vars__.emplace_back(sigma);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, stan::require_std_vector_t<VecVar>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void transform_inits_impl(const stan::io::var_context& context__,
                                   VecI& params_i__, VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.clear();
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      Eigen::Matrix<double, -1, 1> OD_init;
      OD_init = Eigen::Matrix<double, -1, 1>(J);
      stan::math::fill(OD_init, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> OD_init_flat__;
        current_statement__ = 1;
        assign(OD_init_flat__, nil_index_list(), context__.vals_r("OD_init"),
          "assigning variable OD_init_flat__");
        current_statement__ = 1;
        pos__ = 1;
        current_statement__ = 1;
        for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
          current_statement__ = 1;
          assign(OD_init, cons_list(index_uni(sym1__), nil_index_list()),
            OD_init_flat__[(pos__ - 1)], "assigning variable OD_init");
          current_statement__ = 1;
          pos__ = (pos__ + 1);}
      }
      Eigen::Matrix<double, -1, 1> OD_init_free__;
      OD_init_free__ = Eigen::Matrix<double, -1, 1>(J);
      stan::math::fill(OD_init_free__, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 1;
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        current_statement__ = 1;
        assign(OD_init_free__,
          cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_free(OD_init[(sym1__ - 1)], 0),
          "assigning variable OD_init_free__");}
      Eigen::Matrix<double, -1, 1> lam;
      lam = Eigen::Matrix<double, -1, 1>(J);
      stan::math::fill(lam, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> lam_flat__;
        current_statement__ = 2;
        assign(lam_flat__, nil_index_list(), context__.vals_r("lam"),
          "assigning variable lam_flat__");
        current_statement__ = 2;
        pos__ = 1;
        current_statement__ = 2;
        for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
          current_statement__ = 2;
          assign(lam, cons_list(index_uni(sym1__), nil_index_list()),
            lam_flat__[(pos__ - 1)], "assigning variable lam");
          current_statement__ = 2;
          pos__ = (pos__ + 1);}
      }
      Eigen::Matrix<double, -1, 1> lam_free__;
      lam_free__ = Eigen::Matrix<double, -1, 1>(J);
      stan::math::fill(lam_free__, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 2;
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        current_statement__ = 2;
        assign(lam_free__, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_free(lam[(sym1__ - 1)], 0),
          "assigning variable lam_free__");}
      double OD_init_mu;
      OD_init_mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      OD_init_mu = context__.vals_r("OD_init_mu")[(1 - 1)];
      double OD_init_mu_free__;
      OD_init_mu_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      OD_init_mu_free__ = stan::math::lb_free(OD_init_mu, 0);
      double OD_init_sigma;
      OD_init_sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      OD_init_sigma = context__.vals_r("OD_init_sigma")[(1 - 1)];
      double OD_init_sigma_free__;
      OD_init_sigma_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      OD_init_sigma_free__ = stan::math::lb_free(OD_init_sigma, 0);
      double lam_mu;
      lam_mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 5;
      lam_mu = context__.vals_r("lam_mu")[(1 - 1)];
      double lam_mu_free__;
      lam_mu_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 5;
      lam_mu_free__ = stan::math::lb_free(lam_mu, 0);
      double lam_sigma;
      lam_sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 6;
      lam_sigma = context__.vals_r("lam_sigma")[(1 - 1)];
      double lam_sigma_free__;
      lam_sigma_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 6;
      lam_sigma_free__ = stan::math::lb_free(lam_sigma, 0);
      double sigma;
      sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 7;
      sigma = context__.vals_r("sigma")[(1 - 1)];
      double sigma_free__;
      sigma_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 7;
      sigma_free__ = stan::math::lb_free(sigma, 0);
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        vars__.emplace_back(OD_init_free__[(sym1__ - 1)]);}
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        vars__.emplace_back(lam_free__[(sym1__ - 1)]);}
      vars__.emplace_back(OD_init_mu_free__);
      vars__.emplace_back(OD_init_sigma_free__);
      vars__.emplace_back(lam_mu_free__);
      vars__.emplace_back(lam_sigma_free__);
      vars__.emplace_back(sigma_free__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__.clear();
    names__.emplace_back("OD_init");
    names__.emplace_back("lam");
    names__.emplace_back("OD_init_mu");
    names__.emplace_back("OD_init_sigma");
    names__.emplace_back("lam_mu");
    names__.emplace_back("lam_sigma");
    names__.emplace_back("sigma");
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.clear();
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(J)});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(J)});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "OD_init" + '.' + std::to_string(sym1__));
      }}
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "lam" + '.' + std::to_string(sym1__));
      }}
    param_names__.emplace_back(std::string() + "OD_init_mu");
    param_names__.emplace_back(std::string() + "OD_init_sigma");
    param_names__.emplace_back(std::string() + "lam_mu");
    param_names__.emplace_back(std::string() + "lam_sigma");
    param_names__.emplace_back(std::string() + "sigma");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "OD_init" + '.' + std::to_string(sym1__));
      }}
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "lam" + '.' + std::to_string(sym1__));
      }}
    param_names__.emplace_back(std::string() + "OD_init_mu");
    param_names__.emplace_back(std::string() + "OD_init_sigma");
    param_names__.emplace_back(std::string() + "lam_mu");
    param_names__.emplace_back(std::string() + "lam_sigma");
    param_names__.emplace_back(std::string() + "sigma");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"OD_init\",\"type\":{\"name\":\"vector\",\"length\":" << J << "},\"block\":\"parameters\"},{\"name\":\"lam\",\"type\":{\"name\":\"vector\",\"length\":" << J << "},\"block\":\"parameters\"},{\"name\":\"OD_init_mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"OD_init_sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"lam_mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"lam_sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"OD_init\",\"type\":{\"name\":\"vector\",\"length\":" << J << "},\"block\":\"parameters\"},{\"name\":\"lam\",\"type\":{\"name\":\"vector\",\"length\":" << J << "},\"block\":\"parameters\"},{\"name\":\"OD_init_mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"OD_init_sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"lam_mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"lam_sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      std::vector<double> vars_vec(vars.size());
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i) {
        vars.coeffRef(i) = vars_vec[i];
      }
    }

    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }

    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }

    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
  

    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits_impl(context, params_i, params_r_vec, pstream);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i) {
        params_r.coeffRef(i) = params_r_vec[i];
      }
    }
    inline void transform_inits(const stan::io::var_context& context,
                                std::vector<int>& params_i,
                                std::vector<double>& vars,
                                std::ostream* pstream = nullptr) const final {
      transform_inits_impl(context, params_i, vars, pstream);
    }        

};
}

using stan_model = hierarchical_growth_rate_model_namespace::hierarchical_growth_rate_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

stan::math::profile_map& get_stan_profile_data() {
  return hierarchical_growth_rate_model_namespace::profiles__;
}

#endif


