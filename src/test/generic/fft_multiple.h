
/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


//--


#ifndef TEST_GENERIC_FFT_MULTIPLE_H_
#define TEST_GENERIC_FFT_MULTIPLE_H_

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/generic/fft_interface_f.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"


namespace test {


class IOProcFixture : private boost::noncopyable {
 public:
  static const eckit::Configuration & test() {return *getInstance().test_;}

 private:
  static IOProcFixture & getInstance() {
    static IOProcFixture fixture_io_proc;
    return fixture_io_proc;
  }

  IOProcFixture() {
    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "test_fft"));
  }

  ~IOProcFixture() {}

  std::unique_ptr<const eckit::LocalConfiguration> test_;
};


// test #01: validate the input data; checks on data format
void test01() {
  typedef IOProcFixture Test_;

  //--

  const eckit::Configuration * config = &Test_::test();

  // number of sequences
  const int no_seq(config->getInt("no_seq"));
  // number of elements per sequence (excluding the zeros)
  const int no_el_seq(config->getInt("no_el_seq"));
  // total number of elements (including the zeros)
  const int no_el(config->getInt("no_el"));
  // data structure containing the data to be processed/tranformed
  std::vector<float> ds_mfft(config->getFloatVector("ds01"));

  //
  // The following requirements need to be met in order to use the Temperton procedure:
  // 1) The length of the sequence (no_el_seq) must be an even number that has no other
  //    factors except possibly power of 2, 3, and 5.
  // 2) The data structure to be used for the calculations must contain at least 2
  //    sequences (no_seq > 1); all the sequences must have the same number of elements;
  //    each sequence must be followed by 2 extra zeros.
  //    Thus, the size of the data structure must be no_el = (no_seq x (no_el_seq + 2)).
  //

  EXPECT(no_seq > 1);
  EXPECT(no_el_seq % 2 == 0);
  EXPECT((no_seq * (no_el_seq + 2)) == no_el);
  EXPECT(static_cast<int>(ds_mfft.size()) > (2 * (2 + 2)));
  EXPECT(static_cast<int>(ds_mfft.size()) == no_el);

  // number of extra zeros
  const int no_zeros = no_seq*2;

  int idx_z = 0;
  int i = 0;
  while (i < no_zeros) {
    idx_z = idx_z + no_el_seq;
    EXPECT(ds_mfft[idx_z] == 0.0);
    idx_z = idx_z + 1;
    ++i;
    EXPECT(ds_mfft[idx_z] == 0.0);
    idx_z = idx_z + 1;
    ++i;
  }
}

// test #02: evaluate the forward and inverse transforms sequentially
void test02() {
  typedef IOProcFixture Test_;

  //--

  const eckit::Configuration * config = &Test_::test();

  const int no_seq(config->getInt("no_seq"));
  const int no_el_seq(config->getInt("no_el_seq"));
  const int no_el(config->getInt("no_el"));
  std::vector<float> ds_mfft(config->getFloatVector("ds01"));

  std::vector<float> ds_mfft_tmp = ds_mfft;

  // performing forward FFTs, from gridpoint to spectral
  fft_gp2spe(&ds_mfft[0], &no_el, &no_seq, &no_el_seq);

  // performing inverse FFTs, from spectral to gridpoint
  fft_spe2gp(&ds_mfft[0], &no_el, &no_seq, &no_el_seq);

  // tolerance
  const float tol = 1.0e-6;

  // checking that the original sequences have been restored
  for (int i=0; i < no_el; ++i) {
    EXPECT(oops::is_close_absolute(ds_mfft.at(i), ds_mfft_tmp.at(i), tol));
  }
}

// test #03: validate the forward transforms against reference data
void test03() {
  typedef IOProcFixture Test_;

  //--

  const eckit::Configuration * config = &Test_::test();

  const int no_seq(config->getInt("no_seq"));
  const int no_el_seq(config->getInt("no_el_seq"));
  const int no_el(config->getInt("no_el"));
  std::vector<float> ds_mfft(config->getFloatVector("ds01"));

  // performing multiple FFTs (from gridpoint to spectral)
  fft_gp2spe(&ds_mfft[0], &no_el, &no_seq, &no_el_seq);

  // comparing the forward transforms to pre-calculated values
  // (reference data stored in file 'fft_multiple.test')
  for (int i=0; i < no_el; ++i) {
    oops::Log::test() << ds_mfft.at(i) << " ";
  }
  oops::Log::test() << std::flush;
}

// test #04: evaluate multiple FFTs of a set of discrete delta functions
//           (Kronecker delta functions); note that the Fourier Transform
//           of a delta function is a constant function;
void test04() {
  //
  // to prepare the input data ...
  //  1) define 2 discrete delta functions, each function is defined
  //     by a sequence of 6 real values;
  //  2) add 2 extra zeros at the end of each sequence;
  //  3) store the 16 real values into a data structure;
  //

  // number of sequences
  const int no_seq = 2;
  // number of elements per sequence (excluding the zeros)
  const int no_el_seq = 6;
  // total number of elements (including the zeros)
  const int no_el = 16;

  std::vector<float> ds_mfft(no_el, 0.0);
  ds_mfft[0] = 1;
  ds_mfft[8] = 1;

  // performing multiple FFTs (from gridpoint to spectral)
  fft_gp2spe(&ds_mfft[0], &no_el, &no_seq, &no_el_seq);

  // checking that all the real parts of the FFTs coincide and
  // all imaginary parts are equal to zero - i.e. checking that
  // the FFTs are constant functions as expected
  float fft_ref = ds_mfft[0];
  int i = 0;
  while (i < no_el) {
    EXPECT(ds_mfft[i] == fft_ref);
    ++i;
    EXPECT(ds_mfft[i] == 0.0);
    ++i;
  }

  //
  // note that the test has been designed by taking into account
  // the storage scheme used by the Temperton procedure;
  // the forward transforms stored into 'ds_mfft' are complex sequences
  // in Hermitian form; a sequence in Hermitian form of N complex data
  // values can be represented by only N, rather than 2N, independent
  // real values;
  // for details on the storage scheme used by the Temperton procedure
  // see file 'fft_gpoint2spectral_f.F90';
  //
}


class FFTTestsBatch : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~FFTTestsBatch() {}

 private:
  std::string testid() const override {return "test::FFTTestsBatch";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    // test #01
    ts.emplace_back(CASE("FFTTest01")
      { test01(); });
    // test #02
    ts.emplace_back(CASE("FFTTest02")
      { test02(); });
    // test #03
    ts.emplace_back(CASE("FFTTest03")
      { test03(); });
    // test #04
    ts.emplace_back(CASE("FFTTest04")
      { test04(); });
  }

  void clear() const override {}
};

}  // namespace test


#endif  // TEST_GENERIC_FFT_MULTIPLE_H_
