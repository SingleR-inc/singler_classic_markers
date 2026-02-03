<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.12.0">
  <compound kind="file">
    <name>choose.hpp</name>
    <path>singler_classic_markers/</path>
    <filename>choose_8hpp.html</filename>
    <class kind="struct">singler_classic_markers::ChooseOptions</class>
    <namespace>singler_classic_markers</namespace>
  </compound>
  <compound kind="file">
    <name>singler_classic_markers.hpp</name>
    <path>singler_classic_markers/</path>
    <filename>singler__classic__markers_8hpp.html</filename>
    <namespace>singler_classic_markers</namespace>
  </compound>
  <compound kind="struct">
    <name>singler_classic_markers::ChooseBlockedOptions</name>
    <filename>structsingler__classic__markers_1_1ChooseBlockedOptions.html</filename>
    <member kind="variable">
      <type>std::optional&lt; std::size_t &gt;</type>
      <name>number</name>
      <anchorfile>structsingler__classic__markers_1_1ChooseBlockedOptions.html</anchorfile>
      <anchor>a3bde58739f10e836dcdaf87d3a24070a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>keep_ties</name>
      <anchorfile>structsingler__classic__markers_1_1ChooseBlockedOptions.html</anchorfile>
      <anchor>ad79de5bfcb77f7afa0de3d1bd7d47a11</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>use_minimum</name>
      <anchorfile>structsingler__classic__markers_1_1ChooseBlockedOptions.html</anchorfile>
      <anchor>a5646d0a077c38d2b8caad9e4faef0063</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structsingler__classic__markers_1_1ChooseBlockedOptions.html</anchorfile>
      <anchor>a08873af147b688269ac0fc2b97493bb7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>singler_classic_markers::ChooseOptions</name>
    <filename>structsingler__classic__markers_1_1ChooseOptions.html</filename>
    <member kind="variable">
      <type>std::optional&lt; std::size_t &gt;</type>
      <name>number</name>
      <anchorfile>structsingler__classic__markers_1_1ChooseOptions.html</anchorfile>
      <anchor>a278f6b6d5d01a69eb6db7e1b40cc7256</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>keep_ties</name>
      <anchorfile>structsingler__classic__markers_1_1ChooseOptions.html</anchorfile>
      <anchor>adba68b0fee4c7b0f29afb155cab9292f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structsingler__classic__markers_1_1ChooseOptions.html</anchorfile>
      <anchor>a0c59e8aabcaf3bb54f4667f0fdffc946</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>singler_classic_markers</name>
    <filename>namespacesingler__classic__markers.html</filename>
    <class kind="struct">singler_classic_markers::ChooseBlockedOptions</class>
    <class kind="struct">singler_classic_markers::ChooseOptions</class>
    <member kind="function">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::pair&lt; Index_, Stat_ &gt; &gt; &gt; &gt;</type>
      <name>choose</name>
      <anchorfile>namespacesingler__classic__markers.html</anchorfile>
      <anchor>a29f3016c7e7d31e152a4419e8b5edc32</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Label_ *label, const ChooseOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; Index_ &gt; &gt; &gt;</type>
      <name>choose_index</name>
      <anchorfile>namespacesingler__classic__markers.html</anchorfile>
      <anchor>a02dd6670d70dbc7a563785a95282d93b</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Label_ *label, const ChooseOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::pair&lt; Index_, Stat_ &gt; &gt; &gt; &gt;</type>
      <name>choose_blocked</name>
      <anchorfile>namespacesingler__classic__markers.html</anchorfile>
      <anchor>a5ea0e04676bd8d57ada385292c48854f</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Label_ *label, const Block_ *block, const ChooseBlockedOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; Index_ &gt; &gt; &gt;</type>
      <name>choose_blocked_index</name>
      <anchorfile>namespacesingler__classic__markers.html</anchorfile>
      <anchor>a7b2527ca3c9a462e8ac4640aa00890b8</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;matrix, const Label_ *label, const Block_ *block, const ChooseBlockedOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::size_t</type>
      <name>default_number</name>
      <anchorfile>namespacesingler__classic__markers.html</anchorfile>
      <anchor>a860e14fb031f5e489d1da24c65c02f09</anchor>
      <arglist>(std::size_t num_labels)</arglist>
    </member>
  </compound>
</tagfile>
