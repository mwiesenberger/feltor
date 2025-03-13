<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.12.0" doxygen_gitid="c73f5d30f9e8b1df5ba15a1d064ff2067cbb8267">
  <compound kind="file">
    <name>accumulate.cuh</name>
    <path></path>
    <filename>accumulate_8cuh.html</filename>
    <namespace>dg::exblas</namespace>
    <namespace>dg::exblas::gpu</namespace>
  </compound>
  <compound kind="file">
    <name>accumulate.h</name>
    <path></path>
    <filename>accumulate_8h.html</filename>
    <namespace>dg::exblas</namespace>
    <namespace>dg::exblas::cpu</namespace>
  </compound>
  <compound kind="file">
    <name>exdot_cuda.cuh</name>
    <path></path>
    <filename>exdot__cuda_8cuh.html</filename>
    <includes id="accumulate_8cuh" name="accumulate.cuh" local="yes" import="no" module="no" objc="no">accumulate.cuh</includes>
    <namespace>dg::exblas</namespace>
  </compound>
  <compound kind="file">
    <name>exdot_omp.h</name>
    <path></path>
    <filename>exdot__omp_8h.html</filename>
    <includes id="accumulate_8h" name="accumulate.h" local="yes" import="no" module="no" objc="no">accumulate.h</includes>
    <namespace>dg::exblas</namespace>
  </compound>
  <compound kind="file">
    <name>exdot_serial.h</name>
    <path></path>
    <filename>exdot__serial_8h.html</filename>
    <includes id="accumulate_8h" name="accumulate.h" local="yes" import="no" module="no" objc="no">accumulate.h</includes>
    <class kind="union">dg::exblas::udouble</class>
    <class kind="union">dg::exblas::ufloat</class>
    <namespace>dg::exblas</namespace>
  </compound>
  <compound kind="file">
    <name>fpedot_cuda.cuh</name>
    <path></path>
    <filename>fpedot__cuda_8cuh.html</filename>
    <includes id="accumulate_8cuh" name="accumulate.cuh" local="yes" import="no" module="no" objc="no">accumulate.cuh</includes>
    <namespace>dg::exblas</namespace>
  </compound>
  <compound kind="file">
    <name>fpedot_omp.h</name>
    <path></path>
    <filename>fpedot__omp_8h.html</filename>
    <includes id="accumulate_8h" name="accumulate.h" local="yes" import="no" module="no" objc="no">accumulate.h</includes>
    <namespace>dg::exblas</namespace>
  </compound>
  <compound kind="file">
    <name>fpedot_serial.h</name>
    <path></path>
    <filename>fpedot__serial_8h.html</filename>
    <includes id="accumulate_8h" name="accumulate.h" local="yes" import="no" module="no" objc="no">accumulate.h</includes>
    <namespace>dg::exblas</namespace>
  </compound>
  <compound kind="file">
    <name>mpi_accumulate.h</name>
    <path></path>
    <filename>mpi__accumulate_8h.html</filename>
    <includes id="accumulate_8h" name="accumulate.h" local="yes" import="no" module="no" objc="no">accumulate.h</includes>
    <namespace>dg::exblas</namespace>
  </compound>
  <compound kind="union">
    <name>dg::exblas::udouble</name>
    <filename>uniondg_1_1exblas_1_1udouble.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>d</name>
      <anchorfile>uniondg_1_1exblas_1_1udouble.html</anchorfile>
      <anchor>ad18ae17400d59bd70c2de1ed9dd8db43</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int64_t</type>
      <name>i</name>
      <anchorfile>uniondg_1_1exblas_1_1udouble.html</anchorfile>
      <anchor>a14a7a3c5dd02aed2a272ccef78391632</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="union">
    <name>dg::exblas::ufloat</name>
    <filename>uniondg_1_1exblas_1_1ufloat.html</filename>
    <member kind="variable">
      <type>float</type>
      <name>f</name>
      <anchorfile>uniondg_1_1exblas_1_1ufloat.html</anchorfile>
      <anchor>a1a7a1db1feff977d9b2bbe2b161bd070</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int32_t</type>
      <name>i</name>
      <anchorfile>uniondg_1_1exblas_1_1ufloat.html</anchorfile>
      <anchor>af35a6963355371a274b5ac91b6b0599d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>dg::exblas</name>
    <filename>namespacedg_1_1exblas.html</filename>
    <namespace>dg::exblas::cpu</namespace>
    <namespace>dg::exblas::gpu</namespace>
    <class kind="union">dg::exblas::udouble</class>
    <class kind="union">dg::exblas::ufloat</class>
    <member kind="function">
      <type>__host__ void</type>
      <name>exdot_gpu</name>
      <anchorfile>namespacedg_1_1exblas.html</anchorfile>
      <anchor>a48b3fa9c609028c5a26e0280b35fbaa6</anchor>
      <arglist>(unsigned size, PointerOrValue1 x1_ptr, PointerOrValue2 x2_ptr, int64_t *d_superacc, int *status)</arglist>
    </member>
    <member kind="function">
      <type>__host__ void</type>
      <name>exdot_gpu</name>
      <anchorfile>namespacedg_1_1exblas.html</anchorfile>
      <anchor>a95d1caab696c6453b99cdf9ab52daec8</anchor>
      <arglist>(unsigned size, PointerOrValue1 x1_ptr, PointerOrValue2 x2_ptr, PointerOrValue3 x3_ptr, int64_t *d_superacc, int *status)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>exdot_omp</name>
      <anchorfile>namespacedg_1_1exblas.html</anchorfile>
      <anchor>a77b3aef855b6280b5a2fae8e89216e09</anchor>
      <arglist>(unsigned size, PointerOrValue1 x1_ptr, PointerOrValue2 x2_ptr, int64_t *h_superacc, int *status)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>exdot_omp</name>
      <anchorfile>namespacedg_1_1exblas.html</anchorfile>
      <anchor>aa6a5971251a01870e26f2c769fef613d</anchor>
      <arglist>(unsigned size, PointerOrValue1 x1_ptr, PointerOrValue2 x2_ptr, PointerOrValue3 x3_ptr, int64_t *h_superacc, int *status)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>exdot_cpu</name>
      <anchorfile>namespacedg_1_1exblas.html</anchorfile>
      <anchor>a020aba11430835ed502b897fe94eab6d</anchor>
      <arglist>(unsigned size, PointerOrValue1 x1_ptr, PointerOrValue2 x2_ptr, int64_t *h_superacc, int *status)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>exdot_cpu</name>
      <anchorfile>namespacedg_1_1exblas.html</anchorfile>
      <anchor>a5511b744a196e5830b9094b371832bbe</anchor>
      <arglist>(unsigned size, PointerOrValue1 x1_ptr, PointerOrValue2 x2_ptr, PointerOrValue3 x3_ptr, int64_t *h_superacc, int *status)</arglist>
    </member>
    <member kind="function">
      <type>__host__ void</type>
      <name>fpedot_gpu</name>
      <anchorfile>namespacedg_1_1exblas.html</anchorfile>
      <anchor>a87f211c09d47719c9213e789c9cedfdd</anchor>
      <arglist>(int *status, unsigned size, T *fpe, Functor f, PointerOrValues ...xs_ptr)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>fpedot_omp</name>
      <anchorfile>namespacedg_1_1exblas.html</anchorfile>
      <anchor>ad36a221f122d652d75a2d40501d0f9e5</anchor>
      <arglist>(int *status, unsigned size, std::array&lt; T, N &gt; &amp;fpe, Functor f, PointerOrValues ...xs_ptr)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>fpedot_cpu</name>
      <anchorfile>namespacedg_1_1exblas.html</anchorfile>
      <anchor>acfeb0837aac8805f90cbbf6501a41493</anchor>
      <arglist>(int *status, unsigned size, std::array&lt; T, N &gt; &amp;fpe, Functor f, PointerOrValues ...xs_ptr)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>mpi_reduce_communicator</name>
      <anchorfile>namespacedg_1_1exblas.html</anchorfile>
      <anchor>ad475f5ccf3ba2655166b2a6cd90f2475</anchor>
      <arglist>(MPI_Comm comm, MPI_Comm *comm_mod, MPI_Comm *comm_mod_reduce)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reduce_mpi_cpu</name>
      <anchorfile>namespacedg_1_1exblas.html</anchorfile>
      <anchor>aec049a58ea3482c2de1085e8a20bbe0e</anchor>
      <arglist>(unsigned num_superacc, int64_t *in, int64_t *out, MPI_Comm comm, MPI_Comm comm_mod, MPI_Comm comm_mod_reduce)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>dg::exblas::cpu</name>
    <filename>namespacedg_1_1exblas_1_1cpu.html</filename>
    <member kind="function">
      <type>void</type>
      <name>Accumulate</name>
      <anchorfile>namespacedg_1_1exblas_1_1cpu.html</anchorfile>
      <anchor>a2fc13439d5d08d081954a9fe3cb1a5c7</anchor>
      <arglist>(int64_t *accumulator, double x)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>Normalize</name>
      <anchorfile>namespacedg_1_1exblas_1_1cpu.html</anchorfile>
      <anchor>aad356016c718048e6f512642f818f759</anchor>
      <arglist>(int64_t *accumulator, int &amp;imin, int &amp;imax)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>Round</name>
      <anchorfile>namespacedg_1_1exblas_1_1cpu.html</anchorfile>
      <anchor>a8f26aaa4fce96a50841e820f2c1e52e2</anchor>
      <arglist>(int64_t *accumulator)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>dg::exblas::gpu</name>
    <filename>namespacedg_1_1exblas_1_1gpu.html</filename>
    <member kind="function">
      <type>__device__ void</type>
      <name>Accumulate</name>
      <anchorfile>namespacedg_1_1exblas_1_1gpu.html</anchorfile>
      <anchor>acd958aaa7e3e154cb7ac94f3103414cd</anchor>
      <arglist>(int64_t *accumulator, double x, int stride=1)</arglist>
    </member>
    <member kind="function">
      <type>__device__ int</type>
      <name>Normalize</name>
      <anchorfile>namespacedg_1_1exblas_1_1gpu.html</anchorfile>
      <anchor>a37b689b86866fcc9ff825d01c48e1d6c</anchor>
      <arglist>(int64_t *accumulator, int &amp;imin, int &amp;imax, int stride=1)</arglist>
    </member>
    <member kind="function">
      <type>__device__ double</type>
      <name>Round</name>
      <anchorfile>namespacedg_1_1exblas_1_1gpu.html</anchorfile>
      <anchor>a97a51433cb2c261fe77b27bea7b91dc7</anchor>
      <arglist>(int64_t *accumulator)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Extension: ExBLAS</title>
    <filename>index.html</filename>
  </compound>
</tagfile>
