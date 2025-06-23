<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.12.0" doxygen_gitid="c73f5d30f9e8b1df5ba15a1d064ff2067cbb8267">
  <compound kind="file">
    <name>exp_runge_kutta.h</name>
    <path></path>
    <filename>exp__runge__kutta_8h.html</filename>
    <includes id="matrixsqrt_8h" name="matrixsqrt.h" local="yes" import="no" module="no" objc="no">matrixsqrt.h</includes>
    <includes id="lanczos_8h" name="lanczos.h" local="yes" import="no" module="no" objc="no">lanczos.h</includes>
    <includes id="tableau_8h" name="tableau.h" local="yes" import="no" module="no" objc="no">tableau.h</includes>
    <class kind="struct">dg::mat::ExponentialStep</class>
    <class kind="struct">dg::mat::ExponentialERKStep</class>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
  </compound>
  <compound kind="file">
    <name>exp_runge_kutta_t.cpp</name>
    <path></path>
    <filename>exp__runge__kutta__t_8cpp.html</filename>
    <includes id="exp__runge__kutta_8h" name="exp_runge_kutta.h" local="yes" import="no" module="no" objc="no">exp_runge_kutta.h</includes>
    <class kind="struct">MatrixFunction</class>
    <member kind="function">
      <type>void</type>
      <name>rhs</name>
      <anchorfile>exp__runge__kutta__t_8cpp.html</anchorfile>
      <anchor>acfa192a028001deacc8e9fd97695d1e0</anchor>
      <arglist>(double t, double, double &amp;yp)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>solution</name>
      <anchorfile>exp__runge__kutta__t_8cpp.html</anchorfile>
      <anchor>a09689499b4cdb0c1817ac221e6c0dc16</anchor>
      <arglist>(double t, double y0)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>exp__runge__kutta__t_8cpp.html</anchorfile>
      <anchor>ae66f6b31b5ad750f1fe042a706a4e3d4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>matrix</name>
      <anchorfile>exp__runge__kutta__t_8cpp.html</anchorfile>
      <anchor>a1fdef7b422c652726a6b5d351d8d1a4d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>functors.h</name>
    <path></path>
    <filename>functors_8h.html</filename>
    <class kind="struct">dg::mat::BESSELI0</class>
    <class kind="struct">dg::mat::BesselJ</class>
    <class kind="struct">dg::mat::LaguerreL</class>
    <class kind="struct">dg::mat::GAMMA0</class>
    <class kind="struct">dg::mat::GyrolagK</class>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
  </compound>
  <compound kind="file">
    <name>functors_t.cpp</name>
    <path></path>
    <filename>functors__t_8cpp.html</filename>
    <includes id="functors_8h" name="functors.h" local="yes" import="no" module="no" objc="no">functors.h</includes>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>functors__t_8cpp.html</anchorfile>
      <anchor>ae66f6b31b5ad750f1fe042a706a4e3d4</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>invtridiag_t.cpp</name>
    <path></path>
    <filename>invtridiag__t_8cpp.html</filename>
    <includes id="tridiaginv_8h" name="tridiaginv.h" local="yes" import="no" module="no" objc="no">tridiaginv.h</includes>
    <member kind="typedef">
      <type>dg::SquareMatrix&lt; double &gt;</type>
      <name>CooMatrix</name>
      <anchorfile>invtridiag__t_8cpp.html</anchorfile>
      <anchor>a6ede227a90ca6302db81e74afb23dab9</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>dg::TriDiagonal&lt; thrust::host_vector&lt; double &gt; &gt;</type>
      <name>DiaMatrix</name>
      <anchorfile>invtridiag__t_8cpp.html</anchorfile>
      <anchor>a2885407346d145f3c1d0093227a23783</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>dg::HVec</type>
      <name>Container</name>
      <anchorfile>invtridiag__t_8cpp.html</anchorfile>
      <anchor>a6cfea605b3acf549ab5b0948819f83ae</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>invtridiag__t_8cpp.html</anchorfile>
      <anchor>ae66f6b31b5ad750f1fe042a706a4e3d4</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lanczos.h</name>
    <path></path>
    <filename>lanczos_8h.html</filename>
    <includes id="tridiaginv_8h" name="tridiaginv.h" local="yes" import="no" module="no" objc="no">tridiaginv.h</includes>
    <includes id="matrixfunction_8h" name="matrixfunction.h" local="yes" import="no" module="no" objc="no">matrixfunction.h</includes>
    <class kind="class">dg::mat::UniversalLanczos</class>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
  </compound>
  <compound kind="file">
    <name>lanczos_b.cpp</name>
    <path></path>
    <filename>lanczos__b_8cpp.html</filename>
    <includes id="lanczos_8h" name="lanczos.h" local="yes" import="no" module="no" objc="no">lanczos.h</includes>
    <includes id="mcg_8h" name="mcg.h" local="yes" import="no" module="no" objc="no">mcg.h</includes>
    <member kind="typedef">
      <type>dg::DMatrix</type>
      <name>Matrix</name>
      <anchorfile>lanczos__b_8cpp.html</anchorfile>
      <anchor>a319e07af874aea3f863c69aee943e8e4</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>dg::DVec</type>
      <name>Container</name>
      <anchorfile>lanczos__b_8cpp.html</anchorfile>
      <anchor>aea4dfe699d7ea3c32b26265a3cc92b5c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>lhs</name>
      <anchorfile>lanczos__b_8cpp.html</anchorfile>
      <anchor>a74d3e489d871726c04af786a72ca0fdc</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rhs</name>
      <anchorfile>lanczos__b_8cpp.html</anchorfile>
      <anchor>a034dc5c52a5923d4d991998a47b92fbd</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>lanczos__b_8cpp.html</anchorfile>
      <anchor>ae66f6b31b5ad750f1fe042a706a4e3d4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>lx</name>
      <anchorfile>lanczos__b_8cpp.html</anchorfile>
      <anchor>ac058f9d1ca24439eb286a1f7005e5848</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>ly</name>
      <anchorfile>lanczos__b_8cpp.html</anchorfile>
      <anchor>aaa1e0123e1baad493ec687904a428cbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcx</name>
      <anchorfile>lanczos__b_8cpp.html</anchorfile>
      <anchor>ad558bf6ed28fe209c535b30a0b778812</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcy</name>
      <anchorfile>lanczos__b_8cpp.html</anchorfile>
      <anchor>a2581340f1c6a2fcb83593d4fb46d28fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>alpha</name>
      <anchorfile>lanczos__b_8cpp.html</anchorfile>
      <anchor>adeedc97114c5a704906c8e26bdbd2451</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>m</name>
      <anchorfile>lanczos__b_8cpp.html</anchorfile>
      <anchor>a9990e99e87d163c58817550b21d35a83</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>n</name>
      <anchorfile>lanczos__b_8cpp.html</anchorfile>
      <anchor>afd53743d1af6bfa2856255b23b5a41a9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lanczos_mpib.cpp</name>
    <path></path>
    <filename>lanczos__mpib_8cpp.html</filename>
    <includes id="lanczos_8h" name="lanczos.h" local="yes" import="no" module="no" objc="no">lanczos.h</includes>
    <includes id="mcg_8h" name="mcg.h" local="yes" import="no" module="no" objc="no">mcg.h</includes>
    <member kind="typedef">
      <type>dg::MDMatrix</type>
      <name>Matrix</name>
      <anchorfile>lanczos__mpib_8cpp.html</anchorfile>
      <anchor>a541481989f2383d86daa538656d97ff6</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>dg::MDVec</type>
      <name>Container</name>
      <anchorfile>lanczos__mpib_8cpp.html</anchorfile>
      <anchor>a388290df91e467b0a7619c8530016639</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>lhs</name>
      <anchorfile>lanczos__mpib_8cpp.html</anchorfile>
      <anchor>a74d3e489d871726c04af786a72ca0fdc</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rhs</name>
      <anchorfile>lanczos__mpib_8cpp.html</anchorfile>
      <anchor>a034dc5c52a5923d4d991998a47b92fbd</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>lanczos__mpib_8cpp.html</anchorfile>
      <anchor>a0ddf1224851353fc92bfbff6f499fa97</anchor>
      <arglist>(int argc, char *argv[])</arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>lx</name>
      <anchorfile>lanczos__mpib_8cpp.html</anchorfile>
      <anchor>ac058f9d1ca24439eb286a1f7005e5848</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>ly</name>
      <anchorfile>lanczos__mpib_8cpp.html</anchorfile>
      <anchor>aaa1e0123e1baad493ec687904a428cbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcx</name>
      <anchorfile>lanczos__mpib_8cpp.html</anchorfile>
      <anchor>ad558bf6ed28fe209c535b30a0b778812</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcy</name>
      <anchorfile>lanczos__mpib_8cpp.html</anchorfile>
      <anchor>a2581340f1c6a2fcb83593d4fb46d28fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>alpha</name>
      <anchorfile>lanczos__mpib_8cpp.html</anchorfile>
      <anchor>adeedc97114c5a704906c8e26bdbd2451</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>m</name>
      <anchorfile>lanczos__mpib_8cpp.html</anchorfile>
      <anchor>a9990e99e87d163c58817550b21d35a83</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>n</name>
      <anchorfile>lanczos__mpib_8cpp.html</anchorfile>
      <anchor>afd53743d1af6bfa2856255b23b5a41a9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>matrix.h</name>
    <path></path>
    <filename>matrix_8h.html</filename>
    <includes id="functors_8h" name="functors.h" local="yes" import="no" module="no" objc="no">functors.h</includes>
    <includes id="lanczos_8h" name="lanczos.h" local="yes" import="no" module="no" objc="no">lanczos.h</includes>
    <includes id="mcg_8h" name="mcg.h" local="yes" import="no" module="no" objc="no">mcg.h</includes>
    <includes id="matrixsqrt_8h" name="matrixsqrt.h" local="yes" import="no" module="no" objc="no">matrixsqrt.h</includes>
    <includes id="matrixfunction_8h" name="matrixfunction.h" local="yes" import="no" module="no" objc="no">matrixfunction.h</includes>
    <includes id="polarization_8h" name="polarization.h" local="yes" import="no" module="no" objc="no">polarization.h</includes>
    <includes id="polarization__init_8h" name="polarization_init.h" local="yes" import="no" module="no" objc="no">polarization_init.h</includes>
    <includes id="sqrt__cauchy_8h" name="sqrt_cauchy.h" local="yes" import="no" module="no" objc="no">sqrt_cauchy.h</includes>
    <includes id="sqrt__ode_8h" name="sqrt_ode.h" local="yes" import="no" module="no" objc="no">sqrt_ode.h</includes>
    <includes id="tensorelliptic_8h" name="tensorelliptic.h" local="yes" import="no" module="no" objc="no">tensorelliptic.h</includes>
    <includes id="tridiaginv_8h" name="tridiaginv.h" local="yes" import="no" module="no" objc="no">tridiaginv.h</includes>
    <includes id="exp__runge__kutta_8h" name="exp_runge_kutta.h" local="yes" import="no" module="no" objc="no">exp_runge_kutta.h</includes>
  </compound>
  <compound kind="file">
    <name>matrix_doc.h</name>
    <path></path>
    <filename>matrix__doc_8h.html</filename>
  </compound>
  <compound kind="file">
    <name>matrixfunction.h</name>
    <path></path>
    <filename>matrixfunction_8h.html</filename>
    <includes id="functors_8h" name="functors.h" local="yes" import="no" module="no" objc="no">functors.h</includes>
    <includes id="sqrt__cauchy_8h" name="sqrt_cauchy.h" local="yes" import="no" module="no" objc="no">sqrt_cauchy.h</includes>
    <includes id="sqrt__ode_8h" name="sqrt_ode.h" local="yes" import="no" module="no" objc="no">sqrt_ode.h</includes>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
  </compound>
  <compound kind="file">
    <name>matrixfunction_b.cpp</name>
    <path></path>
    <filename>matrixfunction__b_8cpp.html</filename>
    <includes id="lanczos_8h" name="lanczos.h" local="yes" import="no" module="no" objc="no">lanczos.h</includes>
    <includes id="mcg_8h" name="mcg.h" local="yes" import="no" module="no" objc="no">mcg.h</includes>
    <includes id="matrixfunction_8h" name="matrixfunction.h" local="yes" import="no" module="no" objc="no">matrixfunction.h</includes>
    <member kind="function">
      <type>double</type>
      <name>lhs</name>
      <anchorfile>matrixfunction__b_8cpp.html</anchorfile>
      <anchor>a74d3e489d871726c04af786a72ca0fdc</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>matrixfunction__b_8cpp.html</anchorfile>
      <anchor>ae66f6b31b5ad750f1fe042a706a4e3d4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>lx</name>
      <anchorfile>matrixfunction__b_8cpp.html</anchorfile>
      <anchor>ac058f9d1ca24439eb286a1f7005e5848</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>ly</name>
      <anchorfile>matrixfunction__b_8cpp.html</anchorfile>
      <anchor>aaa1e0123e1baad493ec687904a428cbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcx</name>
      <anchorfile>matrixfunction__b_8cpp.html</anchorfile>
      <anchor>ad558bf6ed28fe209c535b30a0b778812</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcy</name>
      <anchorfile>matrixfunction__b_8cpp.html</anchorfile>
      <anchor>a2581340f1c6a2fcb83593d4fb46d28fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>m</name>
      <anchorfile>matrixfunction__b_8cpp.html</anchorfile>
      <anchor>a9990e99e87d163c58817550b21d35a83</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>n</name>
      <anchorfile>matrixfunction__b_8cpp.html</anchorfile>
      <anchor>afd53743d1af6bfa2856255b23b5a41a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>alpha</name>
      <anchorfile>matrixfunction__b_8cpp.html</anchorfile>
      <anchor>adeedc97114c5a704906c8e26bdbd2451</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>ell_fac</name>
      <anchorfile>matrixfunction__b_8cpp.html</anchorfile>
      <anchor>a075e28f6a20c191287b3a283fcc5a01d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>matrixsqrt.h</name>
    <path></path>
    <filename>matrixsqrt_8h.html</filename>
    <includes id="lanczos_8h" name="lanczos.h" local="yes" import="no" module="no" objc="no">lanczos.h</includes>
    <class kind="struct">dg::mat::MatrixSqrt</class>
    <class kind="struct">dg::mat::MatrixFunction</class>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
  </compound>
  <compound kind="file">
    <name>matrixsqrt_b.cpp</name>
    <path></path>
    <filename>matrixsqrt__b_8cpp.html</filename>
    <includes id="lanczos_8h" name="lanczos.h" local="yes" import="no" module="no" objc="no">lanczos.h</includes>
    <includes id="matrixfunction_8h" name="matrixfunction.h" local="yes" import="no" module="no" objc="no">matrixfunction.h</includes>
    <member kind="function">
      <type>double</type>
      <name>lhs</name>
      <anchorfile>matrixsqrt__b_8cpp.html</anchorfile>
      <anchor>a74d3e489d871726c04af786a72ca0fdc</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rhsHelmholtz</name>
      <anchorfile>matrixsqrt__b_8cpp.html</anchorfile>
      <anchor>a3e91244a3704a48cb2e488f015f943f6</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rhsHelmholtzsqrt</name>
      <anchorfile>matrixsqrt__b_8cpp.html</anchorfile>
      <anchor>aa839163482b35e67581302397c3f4995</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>matrixsqrt__b_8cpp.html</anchorfile>
      <anchor>ae66f6b31b5ad750f1fe042a706a4e3d4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>lx</name>
      <anchorfile>matrixsqrt__b_8cpp.html</anchorfile>
      <anchor>ac058f9d1ca24439eb286a1f7005e5848</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>ly</name>
      <anchorfile>matrixsqrt__b_8cpp.html</anchorfile>
      <anchor>aaa1e0123e1baad493ec687904a428cbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcx</name>
      <anchorfile>matrixsqrt__b_8cpp.html</anchorfile>
      <anchor>ad558bf6ed28fe209c535b30a0b778812</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcy</name>
      <anchorfile>matrixsqrt__b_8cpp.html</anchorfile>
      <anchor>a2581340f1c6a2fcb83593d4fb46d28fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>alpha</name>
      <anchorfile>matrixsqrt__b_8cpp.html</anchorfile>
      <anchor>adeedc97114c5a704906c8e26bdbd2451</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>m</name>
      <anchorfile>matrixsqrt__b_8cpp.html</anchorfile>
      <anchor>a9990e99e87d163c58817550b21d35a83</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>n</name>
      <anchorfile>matrixsqrt__b_8cpp.html</anchorfile>
      <anchor>afd53743d1af6bfa2856255b23b5a41a9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>mcg.h</name>
    <path></path>
    <filename>mcg_8h.html</filename>
    <includes id="tridiaginv_8h" name="tridiaginv.h" local="yes" import="no" module="no" objc="no">tridiaginv.h</includes>
    <includes id="matrixfunction_8h" name="matrixfunction.h" local="yes" import="no" module="no" objc="no">matrixfunction.h</includes>
    <class kind="class">dg::mat::MCG</class>
    <class kind="class">dg::mat::MCGFuncEigen</class>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
  </compound>
  <compound kind="file">
    <name>polarization.h</name>
    <path></path>
    <filename>polarization_8h.html</filename>
    <includes id="lanczos_8h" name="lanczos.h" local="yes" import="no" module="no" objc="no">lanczos.h</includes>
    <includes id="matrixsqrt_8h" name="matrixsqrt.h" local="yes" import="no" module="no" objc="no">matrixsqrt.h</includes>
    <includes id="matrixfunction_8h" name="matrixfunction.h" local="yes" import="no" module="no" objc="no">matrixfunction.h</includes>
    <includes id="tensorelliptic_8h" name="tensorelliptic.h" local="yes" import="no" module="no" objc="no">tensorelliptic.h</includes>
    <class kind="class">dg::mat::PolCharge</class>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
  </compound>
  <compound kind="file">
    <name>polarization_b.cpp</name>
    <path></path>
    <filename>polarization__b_8cpp.html</filename>
    <includes id="polarization_8h" name="polarization.h" local="yes" import="no" module="no" objc="no">polarization.h</includes>
    <includes id="functors_8h" name="functors.h" local="yes" import="no" module="no" objc="no">functors.h</includes>
    <member kind="function">
      <type>double</type>
      <name>phi_ana_df</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>a1617cd80c3cb8f7124e11a12fe065b4f</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rho_ana_df</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>a3bdecda61995c1e0e8a11b22500ea432</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>chi_ana</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>ab0c6f9a724136d5d80540f269357027e</anchor>
      <arglist>(double, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rho_ana_FF</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>a10c1d19e3e564501ed9aabc48186d5ed</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>phi_ana_FF</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>a96dc2ced932c5853e97afa38a1affffd</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rho_ana_FFO4</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>a4b71adae5dd78ab72f800967067211a4</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>ae66f6b31b5ad750f1fe042a706a4e3d4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>tau</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>a9a42c4d90eb9808c60bb3f1be6a0c1a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>alpha</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>adeedc97114c5a704906c8e26bdbd2451</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>beta</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>a4049344c4d6020a322197047b513967f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>m</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>a9990e99e87d163c58817550b21d35a83</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>n</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>afd53743d1af6bfa2856255b23b5a41a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>lx</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>ac058f9d1ca24439eb286a1f7005e5848</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>ly</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>aaa1e0123e1baad493ec687904a428cbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcx</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>ad558bf6ed28fe209c535b30a0b778812</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcy</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>a2581340f1c6a2fcb83593d4fb46d28fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>amp</name>
      <anchorfile>polarization__b_8cpp.html</anchorfile>
      <anchor>a18169e9f248b43386608618d6076ae24</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>polarization_init.h</name>
    <path></path>
    <filename>polarization__init_8h.html</filename>
    <class kind="class">dg::mat::PolChargeN</class>
    <class kind="struct">dg::TensorTraits&lt; mat::PolChargeN&lt; G, M, V &gt; &gt;</class>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
  </compound>
  <compound kind="file">
    <name>polarization_init_t.cpp</name>
    <path></path>
    <filename>polarization__init__t_8cpp.html</filename>
    <includes id="polarization__init_8h" name="polarization_init.h" local="yes" import="no" module="no" objc="no">polarization_init.h</includes>
    <class kind="struct">rho_ana</class>
    <member kind="typedef">
      <type>dg::DVec</type>
      <name>SubContainer</name>
      <anchorfile>polarization__init__t_8cpp.html</anchorfile>
      <anchor>a90c9bf21474ba3cda6656d0ae4921cdd</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>phi_ana</name>
      <anchorfile>polarization__init__t_8cpp.html</anchorfile>
      <anchor>abc8fabd103ad2c3d1e1b4e1e8d0ada5b</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>dxphi_ana</name>
      <anchorfile>polarization__init__t_8cpp.html</anchorfile>
      <anchor>a159e407ad3fad968c8ac8c742340bcd6</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>dyphi_ana</name>
      <anchorfile>polarization__init__t_8cpp.html</anchorfile>
      <anchor>a5956b767de9302078b9a2bfc590f831d</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>lapphi_ana</name>
      <anchorfile>polarization__init__t_8cpp.html</anchorfile>
      <anchor>a65af9bbd4129086404984079c570db20</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>polarization__init__t_8cpp.html</anchorfile>
      <anchor>ae66f6b31b5ad750f1fe042a706a4e3d4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>lx</name>
      <anchorfile>polarization__init__t_8cpp.html</anchorfile>
      <anchor>ac058f9d1ca24439eb286a1f7005e5848</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>ly</name>
      <anchorfile>polarization__init__t_8cpp.html</anchorfile>
      <anchor>aaa1e0123e1baad493ec687904a428cbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcx</name>
      <anchorfile>polarization__init__t_8cpp.html</anchorfile>
      <anchor>ad558bf6ed28fe209c535b30a0b778812</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcy</name>
      <anchorfile>polarization__init__t_8cpp.html</anchorfile>
      <anchor>a2581340f1c6a2fcb83593d4fb46d28fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>amp</name>
      <anchorfile>polarization__init__t_8cpp.html</anchorfile>
      <anchor>a18169e9f248b43386608618d6076ae24</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>sqrt_cauchy.h</name>
    <path></path>
    <filename>sqrt__cauchy_8h.html</filename>
    <class kind="struct">dg::mat::SqrtCauchyInt</class>
    <class kind="struct">dg::mat::DirectSqrtCauchy</class>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
    <member kind="define">
      <type>#define</type>
      <name>M_PI</name>
      <anchorfile>sqrt__cauchy_8h.html</anchorfile>
      <anchor>ae71449b1cc6e6250b91f539153a7a0d3</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>sqrt_ode.h</name>
    <path></path>
    <filename>sqrt__ode_8h.html</filename>
    <includes id="tridiaginv_8h" name="tridiaginv.h" local="yes" import="no" module="no" objc="no">tridiaginv.h</includes>
    <class kind="struct">dg::mat::InvSqrtODE</class>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
  </compound>
  <compound kind="file">
    <name>tableau.h</name>
    <path></path>
    <filename>tableau_8h.html</filename>
    <includes id="functors_8h" name="functors.h" local="yes" import="no" module="no" objc="no">functors.h</includes>
    <class kind="struct">dg::mat::FunctionalButcherTableau</class>
    <class kind="struct">dg::mat::ConvertsToFunctionalButcherTableau</class>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
  </compound>
  <compound kind="file">
    <name>tensorelliptic.h</name>
    <path></path>
    <filename>tensorelliptic_8h.html</filename>
    <class kind="struct">dg::mat::TensorElliptic</class>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
  </compound>
  <compound kind="file">
    <name>tensorelliptic2d_b.cpp</name>
    <path></path>
    <filename>tensorelliptic2d__b_8cpp.html</filename>
    <includes id="tensorelliptic_8h" name="tensorelliptic.h" local="yes" import="no" module="no" objc="no">tensorelliptic.h</includes>
    <member kind="function">
      <type>double</type>
      <name>initial</name>
      <anchorfile>tensorelliptic2d__b_8cpp.html</anchorfile>
      <anchor>a8476723b52f7482a3dd190cead682425</anchor>
      <arglist>(double, double)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>pol</name>
      <anchorfile>tensorelliptic2d__b_8cpp.html</anchorfile>
      <anchor>ae3f03ba1bfb09247d2d1650ad8da310b</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rhs</name>
      <anchorfile>tensorelliptic2d__b_8cpp.html</anchorfile>
      <anchor>a034dc5c52a5923d4d991998a47b92fbd</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sol</name>
      <anchorfile>tensorelliptic2d__b_8cpp.html</anchorfile>
      <anchor>a1e3a52b9d799b199247542066121e0b5</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>tensorelliptic2d__b_8cpp.html</anchorfile>
      <anchor>ae66f6b31b5ad750f1fe042a706a4e3d4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>lx</name>
      <anchorfile>tensorelliptic2d__b_8cpp.html</anchorfile>
      <anchor>ac058f9d1ca24439eb286a1f7005e5848</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>ly</name>
      <anchorfile>tensorelliptic2d__b_8cpp.html</anchorfile>
      <anchor>aaa1e0123e1baad493ec687904a428cbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcx</name>
      <anchorfile>tensorelliptic2d__b_8cpp.html</anchorfile>
      <anchor>ad558bf6ed28fe209c535b30a0b778812</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>dg::bc</type>
      <name>bcy</name>
      <anchorfile>tensorelliptic2d__b_8cpp.html</anchorfile>
      <anchor>a2581340f1c6a2fcb83593d4fb46d28fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>amp</name>
      <anchorfile>tensorelliptic2d__b_8cpp.html</anchorfile>
      <anchor>a18169e9f248b43386608618d6076ae24</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>tridiaginv.h</name>
    <path></path>
    <filename>tridiaginv_8h.html</filename>
    <includes id="functors_8h" name="functors.h" local="yes" import="no" module="no" objc="no">functors.h</includes>
    <class kind="class">dg::mat::TridiagInvHMGTI</class>
    <class kind="class">dg::mat::TridiagInvDF</class>
    <class kind="class">dg::mat::TridiagInvD</class>
    <namespace>dg</namespace>
    <namespace>dg::mat</namespace>
  </compound>
  <compound kind="file">
    <name>tridiaginv_b.cpp</name>
    <path></path>
    <filename>tridiaginv__b_8cpp.html</filename>
    <includes id="tridiaginv_8h" name="tridiaginv.h" local="yes" import="no" module="no" objc="no">tridiaginv.h</includes>
    <member kind="typedef">
      <type>double</type>
      <name>value_type</name>
      <anchorfile>tridiaginv__b_8cpp.html</anchorfile>
      <anchor>a099c90519d34c2dbab58e4a771aa6f25</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>dg::SquareMatrix&lt; double &gt;</type>
      <name>CooMatrix</name>
      <anchorfile>tridiaginv__b_8cpp.html</anchorfile>
      <anchor>a6ede227a90ca6302db81e74afb23dab9</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>dg::TriDiagonal&lt; thrust::host_vector&lt; double &gt; &gt;</type>
      <name>DiaMatrix</name>
      <anchorfile>tridiaginv__b_8cpp.html</anchorfile>
      <anchor>a2885407346d145f3c1d0093227a23783</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mu</name>
      <anchorfile>tridiaginv__b_8cpp.html</anchorfile>
      <anchor>a76736576fe2811c803c32143c1e2cdfc</anchor>
      <arglist>(double s, unsigned i, unsigned n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>tridiaginv__b_8cpp.html</anchorfile>
      <anchor>ae66f6b31b5ad750f1fe042a706a4e3d4</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::BESSELI0</name>
    <filename>structdg_1_1mat_1_1_b_e_s_s_e_l_i0.html</filename>
    <templarg>class T</templarg>
    <member kind="function">
      <type></type>
      <name>BESSELI0</name>
      <anchorfile>structdg_1_1mat_1_1_b_e_s_s_e_l_i0.html</anchorfile>
      <anchor>a7603389bb467b5094381e2b3aede9827</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structdg_1_1mat_1_1_b_e_s_s_e_l_i0.html</anchorfile>
      <anchor>a114d1d2ac6f5e8aa68478e8e334512a5</anchor>
      <arglist>(T x) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::BesselJ</name>
    <filename>structdg_1_1mat_1_1_bessel_j.html</filename>
    <templarg>class T</templarg>
    <member kind="function">
      <type></type>
      <name>BesselJ</name>
      <anchorfile>structdg_1_1mat_1_1_bessel_j.html</anchorfile>
      <anchor>acb1e5c85d21919d6a08438670961ce6b</anchor>
      <arglist>(unsigned n)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structdg_1_1mat_1_1_bessel_j.html</anchorfile>
      <anchor>ac11d4de46a26990cffdc9caa0c8111a6</anchor>
      <arglist>(T x) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::ConvertsToFunctionalButcherTableau</name>
    <filename>structdg_1_1mat_1_1_converts_to_functional_butcher_tableau.html</filename>
    <templarg>class real_type</templarg>
    <member kind="typedef">
      <type>real_type</type>
      <name>value_type</name>
      <anchorfile>structdg_1_1mat_1_1_converts_to_functional_butcher_tableau.html</anchorfile>
      <anchor>a38ab07a2ccb81e4a82dacbf87b52e5ad</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ConvertsToFunctionalButcherTableau</name>
      <anchorfile>structdg_1_1mat_1_1_converts_to_functional_butcher_tableau.html</anchorfile>
      <anchor>a2e0af646a1b1e18d9d939e2f8517bb8c</anchor>
      <arglist>(FunctionalButcherTableau&lt; real_type &gt; tableau)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ConvertsToFunctionalButcherTableau</name>
      <anchorfile>structdg_1_1mat_1_1_converts_to_functional_butcher_tableau.html</anchorfile>
      <anchor>a2302b8be69622c413f016e1545b7b45e</anchor>
      <arglist>(enum tableau_identifier id)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ConvertsToFunctionalButcherTableau</name>
      <anchorfile>structdg_1_1mat_1_1_converts_to_functional_butcher_tableau.html</anchorfile>
      <anchor>ab940ba76d4e194145e788c80b814d02a</anchor>
      <arglist>(std::string name)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ConvertsToFunctionalButcherTableau</name>
      <anchorfile>structdg_1_1mat_1_1_converts_to_functional_butcher_tableau.html</anchorfile>
      <anchor>a972d89ed777df58f719fdfa1a5978023</anchor>
      <arglist>(const char *name)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>operator FunctionalButcherTableau&lt; real_type &gt;</name>
      <anchorfile>structdg_1_1mat_1_1_converts_to_functional_butcher_tableau.html</anchorfile>
      <anchor>a08f294f2074e37ab744635b7a755b96e</anchor>
      <arglist>() const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::DirectSqrtCauchy</name>
    <filename>structdg_1_1mat_1_1_direct_sqrt_cauchy.html</filename>
    <templarg>class Container</templarg>
    <member kind="typedef">
      <type>Container</type>
      <name>container_type</name>
      <anchorfile>structdg_1_1mat_1_1_direct_sqrt_cauchy.html</anchorfile>
      <anchor>a4ac289ce8abb2b07e3ec796b80364ce4</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>dg::get_value_type&lt; Container &gt;</type>
      <name>value_type</name>
      <anchorfile>structdg_1_1mat_1_1_direct_sqrt_cauchy.html</anchorfile>
      <anchor>aff9447987a456910589ade8adbaf871d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DirectSqrtCauchy</name>
      <anchorfile>structdg_1_1mat_1_1_direct_sqrt_cauchy.html</anchorfile>
      <anchor>acfe27470b5a9116cd558c2f8c8a6ac26</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DirectSqrtCauchy</name>
      <anchorfile>structdg_1_1mat_1_1_direct_sqrt_cauchy.html</anchorfile>
      <anchor>ac0128f6a44ade23414c2b986f1d54d17</anchor>
      <arglist>(MatrixType &amp;A, const Container &amp;weights, value_type epsCG, unsigned iterCauchy, std::array&lt; value_type, 2 &gt; EVs, int exp)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>construct</name>
      <anchorfile>structdg_1_1mat_1_1_direct_sqrt_cauchy.html</anchorfile>
      <anchor>ac2ffdb9ec63d4bb78afa819709edfdde</anchor>
      <arglist>(Params &amp;&amp;...ps)</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>operator()</name>
      <anchorfile>structdg_1_1mat_1_1_direct_sqrt_cauchy.html</anchorfile>
      <anchor>a87c4810ca8a174e1715b5a04da1e8985</anchor>
      <arglist>(const Container &amp;b, Container &amp;x)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::ExponentialERKStep</name>
    <filename>structdg_1_1mat_1_1_exponential_e_r_k_step.html</filename>
    <templarg>class ContainerType</templarg>
    <member kind="typedef">
      <type>get_value_type&lt; ContainerType &gt;</type>
      <name>value_type</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_e_r_k_step.html</anchorfile>
      <anchor>ad30d6b40b6bafcf24c136663f9c08ea1</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>ContainerType</type>
      <name>container_type</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_e_r_k_step.html</anchorfile>
      <anchor>abe5c01c8a997ba10052856d76d76e622</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ExponentialERKStep</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_e_r_k_step.html</anchorfile>
      <anchor>aae90cab0f8ce631c2f6560dae9614345</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ExponentialERKStep</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_e_r_k_step.html</anchorfile>
      <anchor>a7f6a593dacbc2e0eb2ba81ea9512d56f</anchor>
      <arglist>(ConvertsToFunctionalButcherTableau&lt; value_type &gt; tableau, const ContainerType &amp;copyable)</arglist>
    </member>
    <member kind="function">
      <type>const ContainerType &amp;</type>
      <name>copyable</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_e_r_k_step.html</anchorfile>
      <anchor>a57dc6258958de3a579e6bd90c7cc5917</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>step</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_e_r_k_step.html</anchorfile>
      <anchor>a5389132624790a4ffcbf136029a83dbe</anchor>
      <arglist>(const std::tuple&lt; ExplicitRHS, MatrixFunction &gt; &amp;ode, value_type t0, const ContainerType &amp;u0, value_type &amp;t1, ContainerType &amp;u1, value_type dt, ContainerType &amp;delta)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>step</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_e_r_k_step.html</anchorfile>
      <anchor>a7425ddfb6b534c33ccd0778acf766939</anchor>
      <arglist>(const std::tuple&lt; ExplicitRHS, MatrixFunction &gt; &amp;ode, value_type t0, const ContainerType &amp;u0, value_type &amp;t1, ContainerType &amp;u1, value_type dt)</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>order</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_e_r_k_step.html</anchorfile>
      <anchor>ae540413a2e15bef2dd492c99cd60be05</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>embedded_order</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_e_r_k_step.html</anchorfile>
      <anchor>afcba1ac1f0b2969dba2a23137fb6bc01</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>num_stages</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_e_r_k_step.html</anchorfile>
      <anchor>a7b0c7b082d9c334b7ca9c07881e7949e</anchor>
      <arglist>() const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::ExponentialStep</name>
    <filename>structdg_1_1mat_1_1_exponential_step.html</filename>
    <templarg>class ContainerType</templarg>
    <member kind="typedef">
      <type>get_value_type&lt; ContainerType &gt;</type>
      <name>value_type</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_step.html</anchorfile>
      <anchor>a2ee9727057f99ea845e07e99c3f26b8c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>ContainerType</type>
      <name>container_type</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_step.html</anchorfile>
      <anchor>add5305320388889e0895f8196ba1433e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ExponentialStep</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_step.html</anchorfile>
      <anchor>aabf6e66073cfcda32e071b55975c5810</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ExponentialStep</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_step.html</anchorfile>
      <anchor>a2e92774a8074d5f0ca3875e0e32a16eb</anchor>
      <arglist>(const ContainerType &amp;copyable)</arglist>
    </member>
    <member kind="function">
      <type>const ContainerType &amp;</type>
      <name>copyable</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_step.html</anchorfile>
      <anchor>a4c3d5a4825e99aed070b3a9000e5b5b9</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>step</name>
      <anchorfile>structdg_1_1mat_1_1_exponential_step.html</anchorfile>
      <anchor>a173164b7ba342f14c8335e735bf42c32</anchor>
      <arglist>(MatrixFunction &amp;ode, value_type t0, const ContainerType &amp;u0, value_type &amp;t1, ContainerType &amp;u1, value_type dt)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::FunctionalButcherTableau</name>
    <filename>structdg_1_1mat_1_1_functional_butcher_tableau.html</filename>
    <templarg>class real_type</templarg>
    <member kind="typedef">
      <type>real_type</type>
      <name>value_type</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>ad468425b562cdc90ad6f9fc675032409</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>std::function&lt; value_type(value_type)&gt;</type>
      <name>function_type</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>ad8e39bd1dd2863e6c6af37870187d14f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>FunctionalButcherTableau</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>aa3523dcc41d0438980020bc1b295e2e4</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>FunctionalButcherTableau</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>a37c9cf3a9f6cd200a6b7170956a67269</anchor>
      <arglist>(unsigned s, unsigned order, const function_type *a, const function_type *b, const real_type *c)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>FunctionalButcherTableau</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>ab7b060482888345f726fb8cb9e955147</anchor>
      <arglist>(unsigned s, unsigned embedded_order, unsigned order, const function_type *a, const function_type *b, const function_type *bt, const real_type *c)</arglist>
    </member>
    <member kind="function">
      <type>function_type</type>
      <name>a</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>a9366ecf312288f4b20519756e6ea0a20</anchor>
      <arglist>(unsigned i, unsigned j) const</arglist>
    </member>
    <member kind="function">
      <type>real_type</type>
      <name>c</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>a3966b31a350891bff6beab76a19d94c9</anchor>
      <arglist>(unsigned i) const</arglist>
    </member>
    <member kind="function">
      <type>function_type</type>
      <name>b</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>af05478d9ba10d64ea5580a5f4061ef4c</anchor>
      <arglist>(unsigned j) const</arglist>
    </member>
    <member kind="function">
      <type>function_type</type>
      <name>bt</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>a88a8d01a1e52e5cace091b6968ce2602</anchor>
      <arglist>(unsigned j) const</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>num_stages</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>a19a4c881a70fe0a57c747829f963d6af</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>order</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>a2237a06cf91b8ea0c9f6def9717d120e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>embedded_order</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>ae7086ff7d50a1a8db15307e60dc9a31d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isEmbedded</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>ade4f7f4e6b10839304008234ad5f75a4</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isImplicit</name>
      <anchorfile>structdg_1_1mat_1_1_functional_butcher_tableau.html</anchorfile>
      <anchor>a6c7fcb59321f6320d78fa5d9b70e54b7</anchor>
      <arglist>() const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::GAMMA0</name>
    <filename>structdg_1_1mat_1_1_g_a_m_m_a0.html</filename>
    <templarg>class T</templarg>
    <member kind="function">
      <type></type>
      <name>GAMMA0</name>
      <anchorfile>structdg_1_1mat_1_1_g_a_m_m_a0.html</anchorfile>
      <anchor>a380e515f5e611fb1d0c31c885aff39b2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structdg_1_1mat_1_1_g_a_m_m_a0.html</anchorfile>
      <anchor>a707415547db021be3d707b2eb9bece36</anchor>
      <arglist>(T x) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::GyrolagK</name>
    <filename>structdg_1_1mat_1_1_gyrolag_k.html</filename>
    <templarg>class T</templarg>
    <member kind="function">
      <type></type>
      <name>GyrolagK</name>
      <anchorfile>structdg_1_1mat_1_1_gyrolag_k.html</anchorfile>
      <anchor>a4efc5f465d7fa34e84f9099399a1af98</anchor>
      <arglist>(T n, T a)</arglist>
    </member>
    <member kind="function">
      <type>DG_DEVICE T</type>
      <name>operator()</name>
      <anchorfile>structdg_1_1mat_1_1_gyrolag_k.html</anchorfile>
      <anchor>af88b861402925a38844951eb66ff9ea8</anchor>
      <arglist>(T x) const</arglist>
    </member>
    <member kind="function">
      <type>DG_DEVICE T</type>
      <name>operator()</name>
      <anchorfile>structdg_1_1mat_1_1_gyrolag_k.html</anchorfile>
      <anchor>a86f1e186e5f403c2fac8ed4edb815706</anchor>
      <arglist>(T x, T y) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::InvSqrtODE</name>
    <filename>structdg_1_1mat_1_1_inv_sqrt_o_d_e.html</filename>
    <templarg>class Container</templarg>
    <member kind="typedef">
      <type>Container</type>
      <name>container_type</name>
      <anchorfile>structdg_1_1mat_1_1_inv_sqrt_o_d_e.html</anchorfile>
      <anchor>a4f2e5d6d389de1629438c18ab7544846</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>dg::get_value_type&lt; Container &gt;</type>
      <name>value_type</name>
      <anchorfile>structdg_1_1mat_1_1_inv_sqrt_o_d_e.html</anchorfile>
      <anchor>ad0d846ac8757d5ac85c969beac6cff27</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InvSqrtODE</name>
      <anchorfile>structdg_1_1mat_1_1_inv_sqrt_o_d_e.html</anchorfile>
      <anchor>a799c441360dcf20802bc42a16a6970f1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InvSqrtODE</name>
      <anchorfile>structdg_1_1mat_1_1_inv_sqrt_o_d_e.html</anchorfile>
      <anchor>aef4d53ce6d245919d45d94ff2f9a15ce</anchor>
      <arglist>(MatrixType &amp;A, const Container &amp;copyable)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>construct</name>
      <anchorfile>structdg_1_1mat_1_1_inv_sqrt_o_d_e.html</anchorfile>
      <anchor>a0d7b90cb9321e6cdfd8bc56b33eac166</anchor>
      <arglist>(Params &amp;&amp;...ps)</arglist>
    </member>
    <member kind="function">
      <type>const value_type &amp;</type>
      <name>time</name>
      <anchorfile>structdg_1_1mat_1_1_inv_sqrt_o_d_e.html</anchorfile>
      <anchor>af5c020a1ec9aaf19494348ba5a61e762</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>make_operator</name>
      <anchorfile>structdg_1_1mat_1_1_inv_sqrt_o_d_e.html</anchorfile>
      <anchor>af8d7c552ad4de0eb35de009909addafd</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_inverse_operator</name>
      <anchorfile>structdg_1_1mat_1_1_inv_sqrt_o_d_e.html</anchorfile>
      <anchor>afa9cb6bef239b8968f74a9e664760e45</anchor>
      <arglist>(const MatrixType &amp;OpInv)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>structdg_1_1mat_1_1_inv_sqrt_o_d_e.html</anchorfile>
      <anchor>a6d97d9035316b9ce7e2ee51ba860b788</anchor>
      <arglist>(value_type t, const Container &amp;y, Container &amp;yp)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::LaguerreL</name>
    <filename>structdg_1_1mat_1_1_laguerre_l.html</filename>
    <templarg>class T</templarg>
    <member kind="function">
      <type></type>
      <name>LaguerreL</name>
      <anchorfile>structdg_1_1mat_1_1_laguerre_l.html</anchorfile>
      <anchor>ad100756b1e3a2cc4da1d03b9612f5e51</anchor>
      <arglist>(unsigned n)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>operator()</name>
      <anchorfile>structdg_1_1mat_1_1_laguerre_l.html</anchorfile>
      <anchor>a556faff8bb24777175cc590e0ce6a40c</anchor>
      <arglist>(T x) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::MatrixFunction</name>
    <filename>structdg_1_1mat_1_1_matrix_function.html</filename>
    <templarg>class ContainerType</templarg>
    <member kind="typedef">
      <type>ContainerType</type>
      <name>container_type</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_function.html</anchorfile>
      <anchor>a9e750e904f8571d8685d57b7bfeb12ff</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>dg::get_value_type&lt; ContainerType &gt;</type>
      <name>value_type</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_function.html</anchorfile>
      <anchor>a833f95df0c024c355787db92c11762a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MatrixFunction</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_function.html</anchorfile>
      <anchor>aa7866a262aafdf80cc4b09e72d736ae7</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MatrixFunction</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_function.html</anchorfile>
      <anchor>a2308d30b916cbcf3068bafb114bcb4f3</anchor>
      <arglist>(MatrixType &amp;A, const ContainerType &amp;weights, value_type eps_rel, value_type nrmb_correction=1., unsigned max_iter=500, std::function&lt; value_type(value_type)&gt; f_inner=[](value_type x){return x;})</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>construct</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_function.html</anchorfile>
      <anchor>aee5a2c130e4f445e624dd0dfe04b4f5f</anchor>
      <arglist>(Params &amp;&amp;...ps)</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>get_iter</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_function.html</anchorfile>
      <anchor>a6f745e02da9f2a7f9218553b65f8b137</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_benchmark</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_function.html</anchorfile>
      <anchor>abe550851ad37296dc3aaecb1cade0858</anchor>
      <arglist>(bool benchmark, std::string message=&quot;Function&quot;)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_function.html</anchorfile>
      <anchor>a0d9e6f04568d7a01b8a53a3b1e0e1121</anchor>
      <arglist>(UnaryOp f_outer, const ContainerType0 b, ContainerType1 &amp;x)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>MatrixFunction</name>
    <filename>struct_matrix_function.html</filename>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>struct_matrix_function.html</anchorfile>
      <anchor>abb307ff73fcc3af91123c9bb75f96ce9</anchor>
      <arglist>(UnaryOp f, double x, double &amp;y) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::MatrixSqrt</name>
    <filename>structdg_1_1mat_1_1_matrix_sqrt.html</filename>
    <templarg>class ContainerType</templarg>
    <member kind="typedef">
      <type>ContainerType</type>
      <name>container_type</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_sqrt.html</anchorfile>
      <anchor>a50e9bab2c7f9dfc343e202b44c49a6e1</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>dg::get_value_type&lt; ContainerType &gt;</type>
      <name>value_type</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_sqrt.html</anchorfile>
      <anchor>a446da3d9407473e38c1d4291c9a4d349</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MatrixSqrt</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_sqrt.html</anchorfile>
      <anchor>af30f108bfaa595928a82bc0396da62fe</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MatrixSqrt</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_sqrt.html</anchorfile>
      <anchor>aa9a962de19cc80551dc38c90ecc86416</anchor>
      <arglist>(MatrixType &amp;A, int exp, const ContainerType &amp;weights, value_type eps_rel, value_type nrmb_correction=1., unsigned max_iter=500, unsigned cauchy_steps=40)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>construct</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_sqrt.html</anchorfile>
      <anchor>a19ddf836799708303fa198ba8af8611c</anchor>
      <arglist>(Params &amp;&amp;...ps)</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>get_iter</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_sqrt.html</anchorfile>
      <anchor>a46dccefc0d2a771183823baecf42a9ee</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_benchmark</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_sqrt.html</anchorfile>
      <anchor>a5186e095076c5dcfcc27640737af7237</anchor>
      <arglist>(bool benchmark, std::string message=&quot;SQRT&quot;)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>structdg_1_1mat_1_1_matrix_sqrt.html</anchorfile>
      <anchor>a051d87916930be95c3b63aa70f7d50d0</anchor>
      <arglist>(const ContainerType0 b, ContainerType1 &amp;x)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>dg::mat::MCG</name>
    <filename>classdg_1_1mat_1_1_m_c_g.html</filename>
    <templarg>class ContainerType</templarg>
    <member kind="typedef">
      <type>dg::get_value_type&lt; ContainerType &gt;</type>
      <name>value_type</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g.html</anchorfile>
      <anchor>a6f0ca0749c9c796f2def9a36bf83cda2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MCG</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g.html</anchorfile>
      <anchor>aa31f0d1d1f8de49965949af316176dfd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MCG</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g.html</anchorfile>
      <anchor>ad5c17aaf8001837de28c6b8cc71abca6</anchor>
      <arglist>(const ContainerType &amp;copyable, unsigned max_iterations)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_max</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g.html</anchorfile>
      <anchor>a8f168b1c591719ae51feea8fb7588316</anchor>
      <arglist>(unsigned new_max)</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>get_max</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g.html</anchorfile>
      <anchor>a29b95c39f2324441237783d78f109825</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_verbose</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g.html</anchorfile>
      <anchor>a225770cffcf62afba92454392be71ad3</anchor>
      <arglist>(bool verbose)</arglist>
    </member>
    <member kind="function">
      <type>value_type</type>
      <name>get_bnorm</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g.html</anchorfile>
      <anchor>a0fe8e82d5c91e76d478fd0540aadb38f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>construct</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g.html</anchorfile>
      <anchor>a7eb2fbb134fec5468e36b03dc61b23e7</anchor>
      <arglist>(Params &amp;&amp;...ps)</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>get_iter</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g.html</anchorfile>
      <anchor>ae75f1e7860cb4f8ba4e9c5e22be0b2bc</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Ry</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g.html</anchorfile>
      <anchor>a65a6bce8097ba9596705de0ecb765dce</anchor>
      <arglist>(MatrixType &amp;&amp;A, const DiaMatrixType &amp;T, const ContainerType0 &amp;y, ContainerType1 &amp;x, const ContainerType2 &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g.html</anchorfile>
      <anchor>a843b87e5d71fb9ed55e0da8145be2824</anchor>
      <arglist>(MatrixType &amp;&amp;A, const ContainerType0 &amp;b, const ContainerType1 &amp;weights, value_type eps=1e-12, value_type nrmb_correction=1., value_type res_fac=1.)</arglist>
    </member>
    <member kind="function">
      <type>thrust::host_vector&lt; value_type &gt;</type>
      <name>make_e1</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g.html</anchorfile>
      <anchor>a3e310bafc3660e83f4ba0b6aa2e38b68</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>dg::mat::MCGFuncEigen</name>
    <filename>classdg_1_1mat_1_1_m_c_g_func_eigen.html</filename>
    <templarg>class Container</templarg>
    <member kind="typedef">
      <type>dg::get_value_type&lt; Container &gt;</type>
      <name>value_type</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g_func_eigen.html</anchorfile>
      <anchor>a32fbab67cb186ade43806bed0be7711a</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>dg::HVec</type>
      <name>HVec</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g_func_eigen.html</anchorfile>
      <anchor>a72ea3a23b49729ad1687955ff0ab350f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MCGFuncEigen</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g_func_eigen.html</anchorfile>
      <anchor>ad334170c77fc850aef00455f90c421a2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MCGFuncEigen</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g_func_eigen.html</anchorfile>
      <anchor>ae0510d5ce95c8e9e415c6b75cb96465a</anchor>
      <arglist>(const Container &amp;copyable, unsigned max_iterations)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>construct</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g_func_eigen.html</anchorfile>
      <anchor>a9bbd013151ec5cfc762456968d48d3ab</anchor>
      <arglist>(Params &amp;&amp;...ps)</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_m_c_g_func_eigen.html</anchorfile>
      <anchor>afe1a0a2a8293684603bb84ecc8dec1c5</anchor>
      <arglist>(ContainerType0 &amp;x, UnaryOp f, MatrixType &amp;&amp;A, const ContainerType1 &amp;b, const ContainerType2 &amp;weights, value_type eps, value_type nrmb_correction, value_type res_fac)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>dg::mat::PolCharge</name>
    <filename>classdg_1_1mat_1_1_pol_charge.html</filename>
    <templarg>class Geometry</templarg>
    <templarg>class Matrix</templarg>
    <templarg>class Container</templarg>
    <member kind="typedef">
      <type>get_value_type&lt; Container &gt;</type>
      <name>value_type</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>a183e5d497e21262843d6342763b107ad</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PolCharge</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>a01fcf5e1c895d8c1eaee9e1093c93c95</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PolCharge</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>afceb5f6a6a23e6fb8fb16916e77fc6fa</anchor>
      <arglist>(value_type alpha, std::vector&lt; value_type &gt; eps_gamma, const Geometry &amp;g, direction dir=forward, value_type jfactor=1., std::string mode=&quot;df&quot;, bool commute=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PolCharge</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>a007aee63c42fe3b5d60294088760510b</anchor>
      <arglist>(value_type alpha, std::vector&lt; value_type &gt; eps_gamma, const Geometry &amp;g, bc bcx, bc bcy, direction dir=forward, value_type jfactor=1., std::string mode=&quot;df&quot;, bool commute=false)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>construct</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>a8e9ba1fd6b745072b2c04a2152a9c7d0</anchor>
      <arglist>(Params &amp;&amp;...ps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_chi</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>ac86f1f66727cea55a87291db37bec516</anchor>
      <arglist>(const ContainerType0 &amp;sigma)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_chi</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>a55749459fd4e66caa000842981cd2ae7</anchor>
      <arglist>(const SparseTensor&lt; ContainerType0 &gt; &amp;tau)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_iota</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>a8c5c2edbce770f6240654880f58bb8c8</anchor>
      <arglist>(const ContainerType0 &amp;sigma)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_commute</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>af2adc0da5b4afa914db13315c1ebc046</anchor>
      <arglist>(bool commute)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_commute</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>a47be2dc1eb75714fae3d7c13126b2222</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const Container &amp;</type>
      <name>weights</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>a12efe7bb2faeb0311b0a93fe7699f7d5</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const Container &amp;</type>
      <name>precond</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>acf116a4ea4c3334e79f79281998c8d70</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>variation</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>a56b8684f8d6f00365758cced498b8a65</anchor>
      <arglist>(const ContainerType0 &amp;phi, ContainerType1 &amp;varphi)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>a7713d90d20888fcd259eb02897f155d2</anchor>
      <arglist>(const ContainerType0 &amp;x, ContainerType1 &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>symv</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>a6958e219152873d3eccf20a6901d3f74</anchor>
      <arglist>(const ContainerType0 &amp;x, ContainerType1 &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>symv</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge.html</anchorfile>
      <anchor>a56f5c0cb3b8efba98c316e130333dc5a</anchor>
      <arglist>(value_type alpha, const ContainerType0 &amp;x, value_type beta, ContainerType1 &amp;y)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>dg::mat::PolChargeN</name>
    <filename>classdg_1_1mat_1_1_pol_charge_n.html</filename>
    <templarg>class Geometry</templarg>
    <templarg>class Matrix</templarg>
    <templarg>class Container</templarg>
    <member kind="typedef">
      <type>Geometry</type>
      <name>geometry_type</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>a16c6139771eeda35d603fccc58135593</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Matrix</type>
      <name>matrix_type</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>a8065c9087b271a98aef7cef1af7307af</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Container</type>
      <name>container_type</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>a4b767a38f330c0ab089a02ff901cd953</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>get_value_type&lt; Container &gt;</type>
      <name>value_type</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>a9d9c9c8c8848fef11427edfd0909fbe3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PolChargeN</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>a9c6a8440258289674c6eb36b8aa3a9bb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PolChargeN</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>af5b275d9823f762e00f476238b2bd10f</anchor>
      <arglist>(const Geometry &amp;g, direction dir=forward, value_type jfactor=1.)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PolChargeN</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>af24ccac301df9b0a09ca94e139e08f18</anchor>
      <arglist>(const Geometry &amp;g, bc bcx, bc bcy, direction dir=forward, value_type jfactor=1.)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>construct</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>ace931d4b38025bafb0cdf5a2570fbc5f</anchor>
      <arglist>(Params &amp;&amp;...ps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_phi</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>a0e66661a90315383bdf188bd97a8ef0f</anchor>
      <arglist>(const ContainerType0 &amp;phi)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_dxphi</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>a6143097a0807ac306a06fbc07469b9ef</anchor>
      <arglist>(const ContainerType0 &amp;dxphi)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_dyphi</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>a32a66c69c41d7679c93e025c832aedf9</anchor>
      <arglist>(const ContainerType0 &amp;dyphi)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_lapphi</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>aca4d69dd237167054073538eba1fb6ef</anchor>
      <arglist>(const ContainerType0 &amp;lapphi)</arglist>
    </member>
    <member kind="function">
      <type>const Container &amp;</type>
      <name>weights</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>a68a92e57774140cc300d69df04579028</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const Container &amp;</type>
      <name>precond</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>abd3ffc30012e9571a09b8544a629660b</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>a392b6ca9b8edcced908d0352733d0ba5</anchor>
      <arglist>(const ContainerType0 &amp;x, ContainerType1 &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>symv</name>
      <anchorfile>classdg_1_1mat_1_1_pol_charge_n.html</anchorfile>
      <anchor>a9e83d1e8cd88e6f098cf8c2086290482</anchor>
      <arglist>(const ContainerType0 &amp;x, ContainerType1 &amp;y)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>rho_ana</name>
    <filename>structrho__ana.html</filename>
    <member kind="function">
      <type></type>
      <name>rho_ana</name>
      <anchorfile>structrho__ana.html</anchorfile>
      <anchor>ae525eeee3dde97d61dd9276cf709fd64</anchor>
      <arglist>(dg::Cauchy cauchy)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>structrho__ana.html</anchorfile>
      <anchor>ac9cfd8c9a2a1c83fcf24fda5c4fd88e5</anchor>
      <arglist>(double x, double y) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::SqrtCauchyInt</name>
    <filename>structdg_1_1mat_1_1_sqrt_cauchy_int.html</filename>
    <templarg>class Container</templarg>
    <member kind="typedef">
      <type>Container</type>
      <name>container_type</name>
      <anchorfile>structdg_1_1mat_1_1_sqrt_cauchy_int.html</anchorfile>
      <anchor>ae5d4f5d01a54c092456bbb0bc5ee5d05</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>dg::get_value_type&lt; Container &gt;</type>
      <name>value_type</name>
      <anchorfile>structdg_1_1mat_1_1_sqrt_cauchy_int.html</anchorfile>
      <anchor>a4a7816a96d1d3c42d5eeaab3f2e30b2d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SqrtCauchyInt</name>
      <anchorfile>structdg_1_1mat_1_1_sqrt_cauchy_int.html</anchorfile>
      <anchor>acd9549d3210e9182d62f3047af9476ac</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SqrtCauchyInt</name>
      <anchorfile>structdg_1_1mat_1_1_sqrt_cauchy_int.html</anchorfile>
      <anchor>ab9e01d2e44b9f507b76bac3cc74e16fc</anchor>
      <arglist>(const Container &amp;copyable)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>construct</name>
      <anchorfile>structdg_1_1mat_1_1_sqrt_cauchy_int.html</anchorfile>
      <anchor>a0a38d44a50b268b4b47406c5f4a2508e</anchor>
      <arglist>(Params &amp;&amp;...ps)</arglist>
    </member>
    <member kind="function">
      <type>const double &amp;</type>
      <name>w</name>
      <anchorfile>structdg_1_1mat_1_1_sqrt_cauchy_int.html</anchorfile>
      <anchor>ad467a3c37f8cfa429628fe335c116f69</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>make_denominator</name>
      <anchorfile>structdg_1_1mat_1_1_sqrt_cauchy_int.html</anchorfile>
      <anchor>a7c10a3988d6861d75b063900aef497ba</anchor>
      <arglist>(MatrixType &amp;A) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>structdg_1_1mat_1_1_sqrt_cauchy_int.html</anchorfile>
      <anchor>a749ca76709c85f184f0f23f77b7e7daa</anchor>
      <arglist>(MatrixType0 &amp;&amp;A, MatrixType1 &amp;&amp;wAinv, const ContainerType0 &amp;b, ContainerType1 &amp;x, std::array&lt; value_type, 2 &gt; EVs, unsigned steps, int exp=+1)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::mat::TensorElliptic</name>
    <filename>structdg_1_1mat_1_1_tensor_elliptic.html</filename>
    <templarg>class Geometry</templarg>
    <templarg>class Matrix</templarg>
    <templarg>class Container</templarg>
    <member kind="typedef">
      <type>Container</type>
      <name>container_type</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>a17e12c0b9a8ac65e3804fa1c420b1466</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Geometry</type>
      <name>geometry_type</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>a44f446ae5356f5f1d5d26ea65479b7d8</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Matrix</type>
      <name>matrix_type</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>af75463ea4df7ed4e5c27e951c574632c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>get_value_type&lt; Container &gt;</type>
      <name>value_type</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>ad5a081c256c0888240071fb49205570a</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TensorElliptic</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>a585dcb0b1bc8759d3a27ec31c7625da5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TensorElliptic</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>a241fe4af1df69ee5ae10d2da2a3da6cd</anchor>
      <arglist>(const Geometry &amp;g, direction dir=dg::centered, value_type jfactor=1.)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TensorElliptic</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>ab83f6580cade86e401a1f8544122e9bd</anchor>
      <arglist>(const Geometry &amp;g, bc bcx, bc bcy, direction dir=dg::centered, value_type jfactor=1.)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>construct</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>a569394b77a421fe6b9b88330dc5960ca</anchor>
      <arglist>(Params &amp;&amp;...ps)</arglist>
    </member>
    <member kind="function">
      <type>const Container &amp;</type>
      <name>weights</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>a16dbc07540b3f45e40779c9eabbff1c6</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const Container &amp;</type>
      <name>precond</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>aa738b70a5f2f517fecad565b049deb2f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_chi</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>ac6c8d27ee186e2db22b6b6cc9460f3dc</anchor>
      <arglist>(const ContainerType0 &amp;chi)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_iota</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>a47dc52fac30be8dac614c3d9d8609554</anchor>
      <arglist>(const ContainerType0 &amp;iota)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>variation</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>ad1a8c8197905f0b92c0c6bdb63a44205</anchor>
      <arglist>(const Container &amp;phi, const value_type &amp;alpha, const Container &amp;chi, Container &amp;varphi)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>ac03f35cddf1d0c60555a5b51b6d61c79</anchor>
      <arglist>(const ContainerType0 &amp;x, ContainerType1 &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>symv</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>ae4b009478277ff3494508a5a6656fe0f</anchor>
      <arglist>(const ContainerType0 &amp;x, ContainerType1 &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>symv</name>
      <anchorfile>structdg_1_1mat_1_1_tensor_elliptic.html</anchorfile>
      <anchor>a0f7173ed982f6b0d6964b7e174fbd9b7</anchor>
      <arglist>(value_type alpha, const ContainerType0 &amp;x, value_type beta, ContainerType1 &amp;y)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dg::TensorTraits&lt; mat::PolChargeN&lt; G, M, V &gt; &gt;</name>
    <filename>structdg_1_1_tensor_traits_3_01mat_1_1_pol_charge_n_3_01_g_00_01_m_00_01_v_01_4_01_4.html</filename>
    <templarg>class G</templarg>
    <templarg>class M</templarg>
    <templarg>class V</templarg>
    <member kind="typedef">
      <type>get_value_type&lt; V &gt;</type>
      <name>value_type</name>
      <anchorfile>structdg_1_1_tensor_traits_3_01mat_1_1_pol_charge_n_3_01_g_00_01_m_00_01_v_01_4_01_4.html</anchorfile>
      <anchor>addbcdcb434e68a936d7099ab8f0da1f8</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>SelfMadeMatrixTag</type>
      <name>tensor_category</name>
      <anchorfile>structdg_1_1_tensor_traits_3_01mat_1_1_pol_charge_n_3_01_g_00_01_m_00_01_v_01_4_01_4.html</anchorfile>
      <anchor>a5a56cac329fcc161ca887c6c64782a8b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>dg::mat::TridiagInvD</name>
    <filename>classdg_1_1mat_1_1_tridiag_inv_d.html</filename>
    <templarg>class real_type</templarg>
    <member kind="typedef">
      <type>real_type</type>
      <name>value_type</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d.html</anchorfile>
      <anchor>a5415bb7ebd7abc920722131d668b0c64</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TridiagInvD</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d.html</anchorfile>
      <anchor>af8b994f287445996a35538680bfb8f37</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TridiagInvD</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d.html</anchorfile>
      <anchor>ada36120cf23fbe4000d6468eea112665</anchor>
      <arglist>(const thrust::host_vector&lt; real_type &gt; &amp;copyable)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TridiagInvD</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d.html</anchorfile>
      <anchor>a691f2b152e6022cf374d600ebffac65f</anchor>
      <arglist>(unsigned size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resize</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d.html</anchorfile>
      <anchor>a20cdf7f5246b25626a1cbf0335d2a6f6</anchor>
      <arglist>(unsigned new_size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d.html</anchorfile>
      <anchor>a137463d3259d79684b726963ae40753b</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; real_type &gt; &gt; &amp;T, dg::SquareMatrix&lt; real_type &gt; &amp;Tinv)</arglist>
    </member>
    <member kind="function">
      <type>dg::SquareMatrix&lt; real_type &gt;</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d.html</anchorfile>
      <anchor>ad99795cd2f063ad952ad3f933f321e87</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; real_type &gt; &gt; &amp;T)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d.html</anchorfile>
      <anchor>ac027f1d3217d3251d79acaa5ca42b9e5</anchor>
      <arglist>(const ContainerType0 &amp;a, const ContainerType1 &amp;b, const ContainerType2 &amp;c, dg::SquareMatrix&lt; real_type &gt; &amp;Tinv)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>dg::mat::TridiagInvDF</name>
    <filename>classdg_1_1mat_1_1_tridiag_inv_d_f.html</filename>
    <templarg>class real_type</templarg>
    <member kind="typedef">
      <type>real_type</type>
      <name>value_type</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d_f.html</anchorfile>
      <anchor>a9d66b31144927bc36828c39d051e4250</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TridiagInvDF</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d_f.html</anchorfile>
      <anchor>a8acd3ec48c1868d3701a610d40a6d46a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TridiagInvDF</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d_f.html</anchorfile>
      <anchor>a811edb6ab390cc97cdc2bfee513b51e3</anchor>
      <arglist>(const thrust::host_vector&lt; real_type &gt; &amp;copyable)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TridiagInvDF</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d_f.html</anchorfile>
      <anchor>a1a92b8da47ae2e61ee18625703abdc78</anchor>
      <arglist>(unsigned size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resize</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d_f.html</anchorfile>
      <anchor>a112124d364668a0ace6830f61893ec44</anchor>
      <arglist>(unsigned new_size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d_f.html</anchorfile>
      <anchor>a4bbf49fe814800260247d96d170c568b</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; real_type &gt; &gt; &amp;T, dg::SquareMatrix&lt; real_type &gt; &amp;Tinv)</arglist>
    </member>
    <member kind="function">
      <type>dg::SquareMatrix&lt; real_type &gt;</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d_f.html</anchorfile>
      <anchor>a704a9154a996bc9db4fbedbc244ffed9</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; real_type &gt; &gt; &amp;T)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_d_f.html</anchorfile>
      <anchor>af5894ea4ebfbb56057a1cf5b2b001b65</anchor>
      <arglist>(const ContainerType0 &amp;a, const ContainerType1 &amp;b, const ContainerType2 &amp;c, dg::SquareMatrix&lt; real_type &gt; &amp;Tinv)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>dg::mat::TridiagInvHMGTI</name>
    <filename>classdg_1_1mat_1_1_tridiag_inv_h_m_g_t_i.html</filename>
    <templarg>class real_type</templarg>
    <member kind="typedef">
      <type>real_type</type>
      <name>value_type</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_h_m_g_t_i.html</anchorfile>
      <anchor>a6b29c40b9b593f79fd11b05dae7fd676</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TridiagInvHMGTI</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_h_m_g_t_i.html</anchorfile>
      <anchor>a72a60a8f1de8f392de8073b2261059af</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TridiagInvHMGTI</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_h_m_g_t_i.html</anchorfile>
      <anchor>a0623d58ebd1e083f0c0450ada8552588</anchor>
      <arglist>(const thrust::host_vector&lt; real_type &gt; &amp;copyable)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TridiagInvHMGTI</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_h_m_g_t_i.html</anchorfile>
      <anchor>a11428afeb3d7c23ac2a94ae705d988c8</anchor>
      <arglist>(unsigned size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resize</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_h_m_g_t_i.html</anchorfile>
      <anchor>a7ea38e8ae189f61969dc6ef370dcd126</anchor>
      <arglist>(unsigned new_size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_h_m_g_t_i.html</anchorfile>
      <anchor>a1c624b69336c3773ed3e18812e20e0d1</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; real_type &gt; &gt; &amp;T, dg::SquareMatrix&lt; real_type &gt; &amp;Tinv)</arglist>
    </member>
    <member kind="function">
      <type>dg::SquareMatrix&lt; real_type &gt;</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_h_m_g_t_i.html</anchorfile>
      <anchor>a6aa6a497c616583be89364331baad89c</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; real_type &gt; &gt; &amp;T)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classdg_1_1mat_1_1_tridiag_inv_h_m_g_t_i.html</anchorfile>
      <anchor>a2f609bd61a9aca490ee7678acfe4036c</anchor>
      <arglist>(const ContainerType0 &amp;a, const ContainerType1 &amp;b, const ContainerType2 &amp;c, dg::SquareMatrix&lt; real_type &gt; &amp;Tinv)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>dg::mat::UniversalLanczos</name>
    <filename>classdg_1_1mat_1_1_universal_lanczos.html</filename>
    <templarg>class ContainerType</templarg>
    <member kind="typedef">
      <type>get_value_type&lt; ContainerType &gt;</type>
      <name>value_type</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>a370305216d927c1cbffcfe8fe58ae685</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>UniversalLanczos</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>aa044680e0ff1a7be7528235965c58363</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>UniversalLanczos</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>a75567bfd54a7f388d3abe111ea25554c</anchor>
      <arglist>(const ContainerType &amp;copyable, unsigned max_iterations)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>construct</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>aeb5ed88c774c43d92c8c9a42cf9a755e</anchor>
      <arglist>(Params &amp;&amp;...ps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_max</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>ad367aeffd08562f5e557b633f51e7339</anchor>
      <arglist>(unsigned new_max)</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>get_max</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>a8c3263e23586f66cb758f8bbc94c3e44</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_verbose</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>a7a4938f6a1f1d7a73852f7ba46a7dba8</anchor>
      <arglist>(bool verbose)</arglist>
    </member>
    <member kind="function">
      <type>value_type</type>
      <name>get_bnorm</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>a2a2911d91ba236e25a10c9e01dff5bc1</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>get_iter</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>a9911ecaf6fd13d1d4b7bcdf1691af335</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned</type>
      <name>solve</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>aac8263710136f1baa907c32179d2ad60</anchor>
      <arglist>(ContainerType0 &amp;x, FuncTe1 f, MatrixType &amp;&amp;A, const ContainerType1 &amp;b, const ContainerType2 &amp;weights, value_type eps, value_type nrmb_correction=1., std::string error_norm=&quot;universal&quot;, value_type res_fac=1., unsigned q=1)</arglist>
    </member>
    <member kind="function">
      <type>const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;</type>
      <name>tridiag</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>a542904b98f6391e9c9d8aeed4d624047</anchor>
      <arglist>(MatrixType &amp;&amp;A, const ContainerType0 &amp;b, const ContainerType1 &amp;weights, value_type eps=1e-12, value_type nrmb_correction=1., std::string error_norm=&quot;universal&quot;, value_type res_fac=1., unsigned q=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>normMbVy</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>a5eec63212c438fbd0f3187a934c7948c</anchor>
      <arglist>(MatrixType &amp;&amp;A, const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;T, const ContainerType0 &amp;y, ContainerType1 &amp;x, const ContainerType2 &amp;b, value_type bnorm)</arglist>
    </member>
    <member kind="function">
      <type>const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;</type>
      <name>tridiag</name>
      <anchorfile>classdg_1_1mat_1_1_universal_lanczos.html</anchorfile>
      <anchor>a361b5dcad0b7e2f96b408be1eaa5499b</anchor>
      <arglist>(UnaryOp f, MatrixType &amp;&amp;A, const ContainerType1 &amp;b, const ContainerType2 &amp;weights, value_type eps, value_type nrmb_correction, std::string error_norm=&quot;residual&quot;, value_type res_fac=1., unsigned q=1)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>dg</name>
    <filename>namespacedg.html</filename>
    <namespace>dg::mat</namespace>
    <class kind="struct">dg::TensorTraits&lt; mat::PolChargeN&lt; G, M, V &gt; &gt;</class>
  </compound>
  <compound kind="namespace">
    <name>dg::mat</name>
    <filename>namespacedg_1_1mat.html</filename>
    <class kind="struct">dg::mat::BESSELI0</class>
    <class kind="struct">dg::mat::BesselJ</class>
    <class kind="struct">dg::mat::ConvertsToFunctionalButcherTableau</class>
    <class kind="struct">dg::mat::DirectSqrtCauchy</class>
    <class kind="struct">dg::mat::ExponentialERKStep</class>
    <class kind="struct">dg::mat::ExponentialStep</class>
    <class kind="struct">dg::mat::FunctionalButcherTableau</class>
    <class kind="struct">dg::mat::GAMMA0</class>
    <class kind="struct">dg::mat::GyrolagK</class>
    <class kind="struct">dg::mat::InvSqrtODE</class>
    <class kind="struct">dg::mat::LaguerreL</class>
    <class kind="struct">dg::mat::MatrixFunction</class>
    <class kind="struct">dg::mat::MatrixSqrt</class>
    <class kind="class">dg::mat::MCG</class>
    <class kind="class">dg::mat::MCGFuncEigen</class>
    <class kind="class">dg::mat::PolCharge</class>
    <class kind="class">dg::mat::PolChargeN</class>
    <class kind="struct">dg::mat::SqrtCauchyInt</class>
    <class kind="struct">dg::mat::TensorElliptic</class>
    <class kind="class">dg::mat::TridiagInvD</class>
    <class kind="class">dg::mat::TridiagInvDF</class>
    <class kind="class">dg::mat::TridiagInvHMGTI</class>
    <class kind="class">dg::mat::UniversalLanczos</class>
    <member kind="enumeration">
      <type></type>
      <name>func_tableau_identifier</name>
      <anchorfile>group__exp__int.html</anchorfile>
      <anchor>ga80cb6891d74fa81113e12a8c732bd273</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>EXPLICIT_EULER_1_1</name>
      <anchorfile>group__exp__int.html</anchorfile>
      <anchor>gga80cb6891d74fa81113e12a8c732bd273a62751045292290b7061b31c71d64335f</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>MIDPOINT_2_2</name>
      <anchorfile>group__exp__int.html</anchorfile>
      <anchor>gga80cb6891d74fa81113e12a8c732bd273a0efbd3de504c862e5ae63362d5b87502</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>CLASSIC_4_4</name>
      <anchorfile>group__exp__int.html</anchorfile>
      <anchor>gga80cb6891d74fa81113e12a8c732bd273a28a1dcd699ba7c83154d628694a41587</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>HOCHBRUCK_3_3_4</name>
      <anchorfile>group__exp__int.html</anchorfile>
      <anchor>gga80cb6891d74fa81113e12a8c732bd273aa04d2aa7af544c4eae1a54a4cd21818c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>phi1</name>
      <anchorfile>namespacedg_1_1mat.html</anchorfile>
      <anchor>a212fab563e07542ecb4a2cae109f35c8</anchor>
      <arglist>(T x)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>phi2</name>
      <anchorfile>namespacedg_1_1mat.html</anchorfile>
      <anchor>aa3e0ccd07afd61bd4d29a32f85fa16ac</anchor>
      <arglist>(T x)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>phi3</name>
      <anchorfile>namespacedg_1_1mat.html</anchorfile>
      <anchor>a3b59f00a78f61b301268240e6f9a75ac</anchor>
      <arglist>(T x)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>phi4</name>
      <anchorfile>namespacedg_1_1mat.html</anchorfile>
      <anchor>a9a70cd1283f9d26eaa6b5e52ed49f33e</anchor>
      <arglist>(T x)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>make_FuncEigen_Te1</name>
      <anchorfile>group__tridiagfunction.html</anchorfile>
      <anchor>ga4a848690246725f758dbfaca0412b881</anchor>
      <arglist>(UnaryOp f)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>make_SqrtCauchy_Te1</name>
      <anchorfile>group__tridiagfunction.html</anchorfile>
      <anchor>ga7a8d18ad1d50cc1ecb8c3cad6048b627</anchor>
      <arglist>(int exp, std::array&lt; value_type, 2 &gt; EVs, unsigned stepsCauchy)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>make_SqrtCauchyEigen_Te1</name>
      <anchorfile>group__tridiagfunction.html</anchorfile>
      <anchor>gab28d4900057b7a9c53db4f80c2c0207d</anchor>
      <arglist>(int exp, std::array&lt; value_type, 2 &gt; EVs, unsigned stepsCauchy)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>make_SqrtODE_Te1</name>
      <anchorfile>group__tridiagfunction.html</anchorfile>
      <anchor>gaa9cefff623540bc46612a9643ece3c63</anchor>
      <arglist>(int exp, std::string tableau, value_type rtol, value_type atol, unsigned &amp;number)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>make_Linear_Te1</name>
      <anchorfile>group__tridiagfunction.html</anchorfile>
      <anchor>ga8309fd03f36c68430aa0492a4cd3a492</anchor>
      <arglist>(int exp)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>make_directODESolve</name>
      <anchorfile>group__matrixfunctionapproximation.html</anchorfile>
      <anchor>gadc96cca132705ee18a2c02dfb1b0cb36</anchor>
      <arglist>(ExplicitRHS &amp;&amp;ode, std::string tableau, value_type epsTimerel, value_type epsTimeabs, unsigned &amp;number, value_type t0=0., value_type t1=1.)</arglist>
    </member>
    <member kind="function">
      <type>InvSqrtODE&lt; Container &gt;</type>
      <name>make_inv_sqrtodeCG</name>
      <anchorfile>group__matrixfunctionapproximation.html</anchorfile>
      <anchor>ga4c9bb92436f2d161ca96ce4205edefe7</anchor>
      <arglist>(Matrix &amp;A, const Preconditioner &amp;P, const Container &amp;weights, dg::get_value_type&lt; Container &gt; epsCG)</arglist>
    </member>
    <member kind="function">
      <type>InvSqrtODE&lt; Container &gt;</type>
      <name>make_inv_sqrtodeTri</name>
      <anchorfile>namespacedg_1_1mat.html</anchorfile>
      <anchor>af203581057d74acbf141557a0f407bc9</anchor>
      <arglist>(const Matrix &amp;TH, const Container &amp;copyable)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>make_expode</name>
      <anchorfile>group__matrixfunctionapproximation.html</anchorfile>
      <anchor>ga01fbaf572f56fac2c92900b7275e1a69</anchor>
      <arglist>(MatrixType &amp;A)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>make_besselI0ode</name>
      <anchorfile>group__matrixfunctionapproximation.html</anchorfile>
      <anchor>ga030adac2335f5353c03af46a7a295875</anchor>
      <arglist>(MatrixType &amp;A)</arglist>
    </member>
    <member kind="function">
      <type>value_type</type>
      <name>compute_Tinv_m1</name>
      <anchorfile>group__matrixinvert.html</anchorfile>
      <anchor>gaa70c7b149816c5994f2e9ddec61cf611</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;T, unsigned size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_Tinv_y</name>
      <anchorfile>group__matrixinvert.html</anchorfile>
      <anchor>ga361a500f499d9cf364f7560b01001096</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;T, thrust::host_vector&lt; value_type &gt; &amp;x, const thrust::host_vector&lt; value_type &gt; &amp;y, value_type a=1., value_type d=0.)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>invert</name>
      <anchorfile>group__matrixinvert.html</anchorfile>
      <anchor>gae2802516185970f99ffc804592b4e597</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;T, dg::SquareMatrix&lt; value_type &gt; &amp;Tinv)</arglist>
    </member>
    <member kind="function">
      <type>dg::SquareMatrix&lt; value_type &gt;</type>
      <name>invert</name>
      <anchorfile>group__matrixinvert.html</anchorfile>
      <anchor>ga80471b850c1cb2e30a03be9ea164685f</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;T)</arglist>
    </member>
    <member kind="function">
      <type>std::array&lt; value_type, 2 &gt;</type>
      <name>compute_extreme_EV</name>
      <anchorfile>group__matrixinvert.html</anchorfile>
      <anchor>ga007766053ae09dc8286d345e30838759</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;T)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>matrixnumerical0</name>
    <title>Level 2: Basic numerical algorithms</title>
    <filename>group__matrixnumerical0.html</filename>
    <subgroup>matrixinvert</subgroup>
    <subgroup>tridiagfunction</subgroup>
    <subgroup>matrixfunctionapproximation</subgroup>
    <subgroup>exp_int</subgroup>
  </compound>
  <compound kind="group">
    <name>matrixinvert</name>
    <title>Inversion of tridiagonal matrices</title>
    <filename>group__matrixinvert.html</filename>
    <class kind="class">dg::mat::TridiagInvHMGTI</class>
    <class kind="class">dg::mat::TridiagInvDF</class>
    <class kind="class">dg::mat::TridiagInvD</class>
    <member kind="function">
      <type>value_type</type>
      <name>dg::mat::compute_Tinv_m1</name>
      <anchorfile>group__matrixinvert.html</anchorfile>
      <anchor>gaa70c7b149816c5994f2e9ddec61cf611</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;T, unsigned size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dg::mat::compute_Tinv_y</name>
      <anchorfile>group__matrixinvert.html</anchorfile>
      <anchor>ga361a500f499d9cf364f7560b01001096</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;T, thrust::host_vector&lt; value_type &gt; &amp;x, const thrust::host_vector&lt; value_type &gt; &amp;y, value_type a=1., value_type d=0.)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dg::mat::invert</name>
      <anchorfile>group__matrixinvert.html</anchorfile>
      <anchor>gae2802516185970f99ffc804592b4e597</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;T, dg::SquareMatrix&lt; value_type &gt; &amp;Tinv)</arglist>
    </member>
    <member kind="function">
      <type>dg::SquareMatrix&lt; value_type &gt;</type>
      <name>dg::mat::invert</name>
      <anchorfile>group__matrixinvert.html</anchorfile>
      <anchor>ga80471b850c1cb2e30a03be9ea164685f</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;T)</arglist>
    </member>
    <member kind="function">
      <type>std::array&lt; value_type, 2 &gt;</type>
      <name>dg::mat::compute_extreme_EV</name>
      <anchorfile>group__matrixinvert.html</anchorfile>
      <anchor>ga007766053ae09dc8286d345e30838759</anchor>
      <arglist>(const dg::TriDiagonal&lt; thrust::host_vector&lt; value_type &gt; &gt; &amp;T)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>tridiagfunction</name>
    <title>Tridiagonal Matrix-functions</title>
    <filename>group__tridiagfunction.html</filename>
    <member kind="function">
      <type>auto</type>
      <name>dg::mat::make_FuncEigen_Te1</name>
      <anchorfile>group__tridiagfunction.html</anchorfile>
      <anchor>ga4a848690246725f758dbfaca0412b881</anchor>
      <arglist>(UnaryOp f)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>dg::mat::make_SqrtCauchy_Te1</name>
      <anchorfile>group__tridiagfunction.html</anchorfile>
      <anchor>ga7a8d18ad1d50cc1ecb8c3cad6048b627</anchor>
      <arglist>(int exp, std::array&lt; value_type, 2 &gt; EVs, unsigned stepsCauchy)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>dg::mat::make_SqrtCauchyEigen_Te1</name>
      <anchorfile>group__tridiagfunction.html</anchorfile>
      <anchor>gab28d4900057b7a9c53db4f80c2c0207d</anchor>
      <arglist>(int exp, std::array&lt; value_type, 2 &gt; EVs, unsigned stepsCauchy)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>dg::mat::make_SqrtODE_Te1</name>
      <anchorfile>group__tridiagfunction.html</anchorfile>
      <anchor>gaa9cefff623540bc46612a9643ece3c63</anchor>
      <arglist>(int exp, std::string tableau, value_type rtol, value_type atol, unsigned &amp;number)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>dg::mat::make_Linear_Te1</name>
      <anchorfile>group__tridiagfunction.html</anchorfile>
      <anchor>ga8309fd03f36c68430aa0492a4cd3a492</anchor>
      <arglist>(int exp)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>matrixfunctionapproximation</name>
    <title>Matrix-functions</title>
    <filename>group__matrixfunctionapproximation.html</filename>
    <class kind="class">dg::mat::UniversalLanczos</class>
    <class kind="struct">dg::mat::MatrixSqrt</class>
    <class kind="struct">dg::mat::MatrixFunction</class>
    <class kind="class">dg::mat::MCG</class>
    <class kind="class">dg::mat::MCGFuncEigen</class>
    <class kind="struct">dg::mat::SqrtCauchyInt</class>
    <class kind="struct">dg::mat::DirectSqrtCauchy</class>
    <class kind="struct">dg::mat::InvSqrtODE</class>
    <member kind="function">
      <type>auto</type>
      <name>dg::mat::make_directODESolve</name>
      <anchorfile>group__matrixfunctionapproximation.html</anchorfile>
      <anchor>gadc96cca132705ee18a2c02dfb1b0cb36</anchor>
      <arglist>(ExplicitRHS &amp;&amp;ode, std::string tableau, value_type epsTimerel, value_type epsTimeabs, unsigned &amp;number, value_type t0=0., value_type t1=1.)</arglist>
    </member>
    <member kind="function">
      <type>InvSqrtODE&lt; Container &gt;</type>
      <name>dg::mat::make_inv_sqrtodeCG</name>
      <anchorfile>group__matrixfunctionapproximation.html</anchorfile>
      <anchor>ga4c9bb92436f2d161ca96ce4205edefe7</anchor>
      <arglist>(Matrix &amp;A, const Preconditioner &amp;P, const Container &amp;weights, dg::get_value_type&lt; Container &gt; epsCG)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>dg::mat::make_expode</name>
      <anchorfile>group__matrixfunctionapproximation.html</anchorfile>
      <anchor>ga01fbaf572f56fac2c92900b7275e1a69</anchor>
      <arglist>(MatrixType &amp;A)</arglist>
    </member>
    <member kind="function">
      <type>auto</type>
      <name>dg::mat::make_besselI0ode</name>
      <anchorfile>group__matrixfunctionapproximation.html</anchorfile>
      <anchor>ga030adac2335f5353c03af46a7a295875</anchor>
      <arglist>(MatrixType &amp;A)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>exp_int</name>
    <title>Exponential integrators</title>
    <filename>group__exp__int.html</filename>
    <class kind="struct">dg::mat::ExponentialStep</class>
    <class kind="struct">dg::mat::ExponentialERKStep</class>
    <class kind="struct">dg::mat::FunctionalButcherTableau</class>
    <class kind="struct">dg::mat::ConvertsToFunctionalButcherTableau</class>
    <member kind="enumeration">
      <type></type>
      <name>dg::mat::func_tableau_identifier</name>
      <anchorfile>group__exp__int.html</anchorfile>
      <anchor>ga80cb6891d74fa81113e12a8c732bd273</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>dg::mat::EXPLICIT_EULER_1_1</name>
      <anchorfile>group__exp__int.html</anchorfile>
      <anchor>gga80cb6891d74fa81113e12a8c732bd273a62751045292290b7061b31c71d64335f</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>dg::mat::MIDPOINT_2_2</name>
      <anchorfile>group__exp__int.html</anchorfile>
      <anchor>gga80cb6891d74fa81113e12a8c732bd273a0efbd3de504c862e5ae63362d5b87502</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>dg::mat::CLASSIC_4_4</name>
      <anchorfile>group__exp__int.html</anchorfile>
      <anchor>gga80cb6891d74fa81113e12a8c732bd273a28a1dcd699ba7c83154d628694a41587</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>dg::mat::HOCHBRUCK_3_3_4</name>
      <anchorfile>group__exp__int.html</anchorfile>
      <anchor>gga80cb6891d74fa81113e12a8c732bd273aa04d2aa7af544c4eae1a54a4cd21818c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>matrixnumerical1</name>
    <title>Level 4: Advanced numerical schemes</title>
    <filename>group__matrixnumerical1.html</filename>
    <subgroup>matrixmatrixoperators</subgroup>
  </compound>
  <compound kind="group">
    <name>matrixmatrixoperators</name>
    <title>Elliptic operators</title>
    <filename>group__matrixmatrixoperators.html</filename>
    <class kind="class">dg::mat::PolCharge</class>
    <class kind="class">dg::mat::PolChargeN</class>
    <class kind="struct">dg::mat::TensorElliptic</class>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Extension: Matrix functions</title>
    <filename>index.html</filename>
  </compound>
</tagfile>
