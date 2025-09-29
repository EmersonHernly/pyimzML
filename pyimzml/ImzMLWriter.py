# TODO: Add a variable to the writer class that allows for saving scan/spectrum groups as referenceable param groups to avoid redundancy.
# TODO: Add support for MSn scans, with multiple precursors. For now, hardcode at one window/selected ion.
# TODO: add support for correct labelling of other scan types such as MRM.


from __future__ import print_function

import os
import numpy as np
import uuid
import hashlib
import sys
import getopt
from collections import namedtuple, OrderedDict, defaultdict

from wheezy.template import Engine, CoreExtension, DictLoader

from pyimzml.compression import NoCompression, ZlibCompression

IMZML_TEMPLATE = """\
@require(uuid, sha1sum, mz_data_type, int_data_type, run_id, spectra, mode, obo_codes, obo_names, mz_compression, int_compression, polarity, spec_type, scan_direction, scan_pattern, scan_type, line_scan_direction, ms_levels, image_x_dimension, image_y_dimension, pixel_size_x, pixel_size_y, xml_element_strings)
<?xml version="1.0" encoding="ISO-8859-1"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd" version="1.1">
  <cvList count="2">
    <cv URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" id="MS" version="3.65.0"/>
    <cv URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo" fullName="Unit Ontology" id="UO" version="12:10:2011"/>
    <cv URI="http://www.maldi-msi.org/download/imzml/imagingMS.obo" fullName="Imaging MS Ontology" id="IMS" version="0.9.1"/>
  </cvList>
  <fileDescription>
    <fileContent>
      @for level in ms_levels:
      @if level==1:
      <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
      @else:
      <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>
      @end
      @end
      @if spec_type=='centroid':
      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
      @elif spec_type=='profile':
      <cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>
      @end
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[mode]" name="@mode" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000080" name="universally unique identifier" value="@uuid"/>
      <cvParam cvRef="IMS" accession="IMS:1000091" name="ibd SHA-1" value="@sha1sum"/>
    </fileContent>
    @if xml_element_strings.get("source_file_list") is not None:
    <sourceFileList count="@{xml_element_strings.get("source_file_count")!!s}">
    @for source_file_string in xml_element_strings.get("source_file_list"):
    @source_file_string
    @end
    </sourceFileList>
    @end
  </fileDescription>
  @if xml_element_strings.get("referenceable_param_group_list_element") is not None:
  <referenceableParamGroupList count="5">
  @for ref_param_string in xml_element_strings["referenceable_param_group_list_element"]:
    @ref_param_string
  @end
  @else:
  <referenceableParamGroupList count="4">
  @end
    <referenceableParamGroup id="mzArray">
      <cvParam cvRef="MS" accession="MS:@obo_codes[mz_compression]" name="@mz_compression" value=""/>
      <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
      <cvParam cvRef="MS" accession="MS:@obo_codes[mz_data_type]" name="@mz_data_type" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>
    </referenceableParamGroup>
    <referenceableParamGroup id="intensityArray">
      <cvParam cvRef="MS" accession="MS:@obo_codes[int_data_type]" name="@int_data_type" value=""/>
      <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>
      <cvParam cvRef="MS" accession="MS:@obo_codes[int_compression]" name="@int_compression" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>
    </referenceableParamGroup>
    <referenceableParamGroup id="scan1">
      <cvParam cvRef="MS" accession="MS:1000093" name="increasing m/z scan"/>
    </referenceableParamGroup>
    <referenceableParamGroup id="spectrum1">
      @if spec_type=='centroid':
      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
      @elif spec_type=='profile':
      <cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>
      @end
      @if polarity=='positive':
      <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>
      @elif polarity=='negative':
      <cvParam cvRef="MS" accession="MS:1000129" name="negative scan" value=""/>
      @end
    </referenceableParamGroup>
  </referenceableParamGroupList>
  @if xml_element_strings.get("software_list_element") is None:
  <softwareList count="1">
  @else:
  <softwareList count="@{xml_element_strings.get("software_list_count")!!s}">
  @for software_string in xml_element_strings["software_list_element"]:
  @software_string
  @end
  @end
    <software id="pyimzml" version="0.0001">
      <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value=""/>
    </software>
  </softwareList>
  <scanSettingsList count="1">
    <scanSettings id="scanSettings1">
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[scan_direction]" name="@obo_names[scan_direction]"/>
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[scan_pattern]" name="@obo_names[scan_pattern]"/>
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[scan_type]" name="@obo_names[scan_type]"/>
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[line_scan_direction]" name="@obo_names[line_scan_direction]"/>
      <cvParam cvRef="IMS" accession="IMS:1000042" name="max count of pixels x" value="@{str(max(s["coords"][0] for s in spectra))!!s}"/>
      <cvParam cvRef="IMS" accession="IMS:1000043" name="max count of pixels y" value="@{str(max(s["coords"][1] for s in spectra))!!s}"/>
      @if image_x_dimension is not None:
      <cvParam cvRef="IMS" accession="IMS:1000044" name="max dimension x" value="@{str(image_x_dimension)!!s}" unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"/>
      <cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size (x)" value="@{str(pixel_size_x)!!s}" unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"/>
      @end
      @if image_y_dimension is not None:
      <cvParam cvRef="IMS" accession="IMS:1000045" name="max dimension y" value="@{str(image_y_dimension)!!s}" unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"/>
      <cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y" value="@{str(pixel_size_y)!!s}" unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"/>
      @end
    </scanSettings>
  </scanSettingsList>
  @if xml_element_strings.get("instrument_configuration_list_element") is None:
  <instrumentConfigurationList count="1">
    <instrumentConfiguration id="IC1">
    </instrumentConfiguration>
  </instrumentConfigurationList>
  @else:
  @for instrument_config in xml_element_strings["instrument_configuration_list_element"]:
  @instrument_config
  @end
  @end
  @if xml_element_strings.get("data_processing_list_element") is None:
  <dataProcessingList count="1">
  @else:
  <dataProcessingList count="@{xml_element_strings.get("data_processing_list_count")!!s}">
  @for data_process in xml_element_strings.get("data_processing_list_element"):
  @data_process
  @end
  @end
    <dataProcessing id="export_from_pyimzml">
      <processingMethod order="0" softwareRef="pyimzml">
        <cvParam cvRef="IMS" accession="IMS:1000500" name="conversion to imzML"/>
      </processingMethod>
    </dataProcessing>
  </dataProcessingList>
  <run defaultInstrumentConfigurationRef="IC1" id="@run_id">
    <spectrumList count="@{len(spectra)!!s}" defaultDataProcessingRef="export_from_pyimzml">
      @for index, s in enumerate(spectra):
      <spectrum defaultArrayLength="0" id="spectrum=@{(index+1)!!s}" index="@{(index+1)!!s}">
        <referenceableParamGroupRef ref="spectrum1"/>
        @if s["ms_level"]==1:
        <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
        @elif s["ms_level"]!=1:
        <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>
        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="@{s["ms_level"]!!s}"/>
        @end
        <cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="@{s["mz_min"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
        <cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="@{s["mz_max"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
        <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="@{s["mz_base"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
        <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="@{s["int_base"]!!s}" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of counts"/>
        <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="@{s["int_tic"]!!s}"/>
        <scanList count="1">
          <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>
          <scan instrumentConfigurationRef="IC1">
            <referenceableParamGroupRef ref="scan1"/>
            @if s.get("scan_start_time") is not None:
            <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="@{s["scan_start_time"]!!s}" unitCvRef="unit.ontology" unitAccession="UO:000031" unitName="minute"/>
            @end
            @if s.get("filter_string"):
            <cvParam cvRef="MS" accession="MS:1000512" name="filter string" value="@{s["filter_string"]!!s}"/>
            @else:
            <cvParam cvRef="MS" accession="MS:1000512" name="filter string" value=""/>
            @end
            <cvParam cvRef="IMS" accession="IMS:1000050" name="position x" value="@{s["coords"][0]!!s}"/>
            <cvParam cvRef="IMS" accession="IMS:1000051" name="position y" value="@{s["coords"][1]!!s}"/>
            @if len(s["coords"]) == 3:
            <cvParam cvRef="IMS" accession="IMS:1000052" name="position z" value="@{s["coords"][2]!!s}"/>
            @end
            @if s["userParams"]:
                @for up in s["userParams"]:
                <userParam name="@up['name']" value="@up['value']"/> 
                @end
            @end
            <scanWindowList count="1">
              <scanWindow>
                <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="@{s["mass_window"][0]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
                <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="@{s["mass_window"][1]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              </scanWindow>
            </scanWindowList>
          </scan>
        </scanList>
        @if s["ms_level"]>1:
        @if s.get("precursor_element_string") is None:
        <precursorList count="1">
          <precursor>
            <isolationWindow>
              <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="@{s["precursor_mz"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              @if s["isolation_window_lower_offset"] is not None:
              <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="@{s["isolation_window_lower_offset"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              @else:
              <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="0.5" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              @end
              @if s["isolation_window_upper_offset"] is not None:
              <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="@{s["isolation_window_upper_offset"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              @else:
              <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.5" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              @end
            </isolationWindow>
            <selectedIonList count="1">
              <selectedIon>
                <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="@{s["precursor_mz"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              </selectedIon>
            </selectedIonList>
            @if s.get("activation") is not None:
            <activation>
              <cvParam cvRef="MS" accession="MS:1000422" name="beam-type collision-induced dissociation" value=""/>
              <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35.0" unitCvRef="unit.ontology" unitAccession="UO:0000266" unitName="electronvolt"/>
            </activation>
            @end
          </precursor>
        </precursorList>
        @else:
        @for precursor_element_line in s["precursor_element_string"]:
        @precursor_element_line
        @end
        @end
        @end
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="mzArray"/>
            <cvParam cvRef="IMS" accession="IMS:1000103" name="external array length" value="@{s["mz_len"]!!s}"/>
            <cvParam cvRef="IMS" accession="IMS:1000104" name="external encoded length" value="@{s["mz_enc_len"]!!s}"/>
            <cvParam cvRef="IMS" accession="IMS:1000102" name="external offset" value="@{s["mz_offset"]!!s}"/>
            <binary/>
          </binaryDataArray>
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="intensityArray"/>
            <cvParam cvRef="IMS" accession="IMS:1000103" name="external array length" value="@{s["int_len"]!!s}"/>
            <cvParam cvRef="IMS" accession="IMS:1000104" name="external encoded length" value="@{s["int_enc_len"]!!s}"/>
            <cvParam cvRef="IMS" accession="IMS:1000102" name="external offset" value="@{s["int_offset"]!!s}"/>
            <binary/>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
      @end
    </spectrumList>
  </run>
</mzML>
"""

IMZML_MOBILITY_TEMPLATE = """\
@require(uuid, sha1sum, mz_data_type, int_data_type, mob_data_type, run_id, spectra, mode, obo_codes, obo_names, mz_compression, int_compression, mob_compression, polarity, spec_type, scan_direction, scan_pattern, scan_type, line_scan_direction, ms_levels, image_x_dimension, image_y_dimension, pixel_size_x, pixel_size_y, mobility_name, mobility_accession, mobility_unit, mobility_unit_accession, xml_element_strings)
<?xml version="1.0" encoding="ISO-8859-1"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd" version="1.1">
  <cvList count="2">
    <cv URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" id="MS" version="3.65.0"/>
    <cv URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo" fullName="Unit Ontology" id="UO" version="12:10:2011"/>
    <cv URI="http://www.maldi-msi.org/download/imzml/imagingMS.obo" fullName="Imaging MS Ontology" id="IMS" version="0.9.1"/>
  </cvList>

  <fileDescription>
    <fileContent>
      @for level in ms_levels:
      @if level==1:
      <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
      @else:
      <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>
      @end
      @end
      @if spec_type=='centroid':
      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
      @elif spec_type=='profile':
      <cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>
      @end
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[mode]" name="@mode" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000080" name="universally unique identifier" value="@uuid"/>
      <cvParam cvRef="IMS" accession="IMS:1000091" name="ibd SHA-1" value="@sha1sum"/>
    </fileContent>
    @if xml_element_strings.get("source_file_list") is not None:
    <sourceFileList count="@{xml_element_strings.get("source_file_count")!!s}">
    @for source_file_string in xml_element_strings.get("source_file_list"):
    @source_file_string
    @end
    </sourceFileList>
    @end

  </fileDescription>

  @if xml_element_strings.get("referenceable_param_group_list_element") is not None:
  <referenceableParamGroupList count="6">
  @for ref_param_string in xml_element_strings["referenceable_param_group_list_element"]:
    @ref_param_string
  @end
  @else:
  <referenceableParamGroupList count="5">
  @end
    <referenceableParamGroup id="mzArray">
      <cvParam cvRef="MS" accession="MS:@obo_codes[mz_compression]" name="@mz_compression" value=""/>
      <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
      <cvParam cvRef="MS" accession="MS:@obo_codes[mz_data_type]" name="@mz_data_type" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>
    </referenceableParamGroup>
    <referenceableParamGroup id="intensityArray">
      <cvParam cvRef="MS" accession="MS:@obo_codes[int_data_type]" name="@int_data_type" value=""/>
      <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>
      <cvParam cvRef="MS" accession="MS:@obo_codes[int_compression]" name="@int_compression" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>
    </referenceableParamGroup>
    <referenceableParamGroup id="mobilityArray">
      <cvParam cvRef="MS" accession="MS:@obo_codes[mob_compression]" name="@mob_compression" value=""/>
      <cvParam cvRef="MS" accession="@{str(mobility_accession)!!s}" name="@{str(mobility_name)!!s}" unitCvRef="MS" unitAccession="@{str(mobility_unit_accession)!!s}" unitName="@{str(mobility_unit)!!s}"/>
      <cvParam cvRef="MS" accession="MS:@obo_codes[mob_data_type]" name="@mob_data_type" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>
    </referenceableParamGroup>
    <referenceableParamGroup id="scan1">
      <cvParam cvRef="MS" accession="MS:1000093" name="increasing m/z scan"/>
    </referenceableParamGroup>
    <referenceableParamGroup id="spectrum1">
      @if spec_type=='centroid':
      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
      @elif spec_type=='profile':
      <cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>
      @end
      @if polarity=='positive':
      <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>
      @elif polarity=='negative':
      <cvParam cvRef="MS" accession="MS:1000129" name="negative scan" value=""/>
      @end
    </referenceableParamGroup>
  </referenceableParamGroupList>

  @if xml_element_strings.get("software_list_element") is not None:
  <softwareList count="1">
  @else:
  <softwareList count="@{xml_element_strings.get("software_list_count")!!s}">
  @for software_string in xml_element_strings["software_list_element"]:
  @software_string
  @end
  @end
    <software id="pyimzml" version="0.0001">
      <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value=""/>
    </software>
  </softwareList>

  <scanSettingsList count="1">
    <scanSettings id="scanSettings1">
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[scan_direction]" name="@obo_names[scan_direction]"/>
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[scan_pattern]" name="@obo_names[scan_pattern]"/>
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[scan_type]" name="@obo_names[scan_type]"/>
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[line_scan_direction]" name="@obo_names[line_scan_direction]"/>
      <cvParam cvRef="IMS" accession="IMS:1000042" name="max count of pixels x" value="@{str(max(s["coords"][0] for s in spectra))!!s}"/>
      <cvParam cvRef="IMS" accession="IMS:1000043" name="max count of pixels y" value="@{str(max(s["coords"][1] for s in spectra))!!s}"/>
      @if image_x_dimension is not None:
      <cvParam cvRef="IMS" accession="IMS:1000044" name="max dimension x" value="@{str(image_x_dimension)!!s}" unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"/>
      <cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size (x)" value="@{str(pixel_size_x)!!s}" unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"/>
      @end
      @if image_y_dimension is not None:
      <cvParam cvRef="IMS" accession="IMS:1000045" name="max dimension y" value="@{str(image_y_dimension)!!s}" unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"/>
      <cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y" value="@{str(pixel_size_y)!!s}" unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"/>
      @end
    </scanSettings>
  </scanSettingsList>

  @if xml_element_strings.get("instrument_configuration_list_element") is None:
  <instrumentConfigurationList count="1">
    <instrumentConfiguration id="IC1">
    </instrumentConfiguration>
  </instrumentConfigurationList>
  @else:
  @for instrument_config in xml_element_strings["instrument_configuration_list_element"]:
  @instrument_config
  @end
  @end
  @if xml_element_strings.get("data_processing_list_element") is None:
  <dataProcessingList count="1">
  @else:
  <dataProcessingList count="@{xml_element_strings.get("data_processing_list_count")!!s}">
  @for data_process in xml_element_strings.get("data_processing_list_element"):
  @data_process
  @end
  @end
    <dataProcessing id="export_from_pyimzml">
      <processingMethod order="0" softwareRef="pyimzml">
        <cvParam cvRef="IMS" accession="IMS:1000500" name="conversion to imzML"/>
      </processingMethod>
    </dataProcessing>
  </dataProcessingList>

  <run defaultInstrumentConfigurationRef="IC1" id="@run_id">
    <spectrumList count="@{len(spectra)!!s}" defaultDataProcessingRef="export_from_pyimzml">
      @for index, s in enumerate(spectra):
      <spectrum defaultArrayLength="0" id="spectrum=@{(index+1)!!s}" index="@{(index+1)!!s}">
        <referenceableParamGroupRef ref="spectrum1"/>
        @if s["ms_level"]==1:
        <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
        @elif s["ms_level"]!=1:
        <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>
        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="@{s["ms_level"]!!s}"/>
        @end
        <cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="@{s["mz_min"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
        <cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="@{s["mz_max"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
        <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="@{s["mz_base"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
        <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="@{s["int_base"]!!s}" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of counts"/>
        <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="@{s["int_tic"]!!s}"/>
        <scanList count="1">
          <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>
          <scan instrumentConfigurationRef="IC1">
            <referenceableParamGroupRef ref="scan1"/>
            @if s.get("scan_start_time") is not None:
            <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="@{s["scan_start_time"]!!s}" unitCvRef="unit.ontology" unitAccession="UO:000031" unitName="minute"/>
            @end
            @if s.get("filter_string"):
            <cvParam cvRef="MS" accession="MS:1000512" name="filter string" value="@{s["filter_string"]!!s}"/>
            @else:
            <cvParam cvRef="MS" accession="MS:1000512" name="filter string" value=""/>
            @end
            <cvParam cvRef="IMS" accession="IMS:1000050" name="position x" value="@{s["coords"][0]!!s}"/>
            <cvParam cvRef="IMS" accession="IMS:1000051" name="position y" value="@{s["coords"][1]!!s}"/>
            @if len(s["coords"]) == 3:
            <cvParam cvRef="IMS" accession="IMS:1000052" name="position z" value="@{s["coords"][2]!!s}"/>
            @end
            @if s["userParams"]:
                @for up in s["userParams"]:
                <userParam name="@up['name']" value="@up['value']"/> 
                @end
            @end
            <scanWindowList count="1">
              <scanWindow>
                <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="@{s["mass_window"][0]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
                <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="@{s["mass_window"][1]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              </scanWindow>
            </scanWindowList>
          </scan>
        </scanList>

        @if s["ms_level"]>1:
        @if s.get("precursor_element_string") is None:
        <precursorList count="1">
          <precursor>
            <isolationWindow>
              <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="@{s["precursor_mz"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              @if s["isolation_window_lower_offset"] is not None:
              <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="@{s["isolation_window_lower_offset"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              @else:
              <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="0.5" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              @end
              @if s["isolation_window_upper_offset"] is not None:
              <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="@{s["isolation_window_upper_offset"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              @else:
              <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.5" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              @end
            </isolationWindow>
            <selectedIonList count="1">
              <selectedIon>
                <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="@{s["precursor_mz"]!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              </selectedIon>
            </selectedIonList>
            @if s.get("activation") is not None:
            <activation>
              <cvParam cvRef="MS" accession="MS:1000422" name="beam-type collision-induced dissociation" value=""/>
              <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35.0" unitCvRef="unit.ontology" unitAccession="UO:0000266" unitName="electronvolt"/>
            </activation>
            @end
          </precursor>
        </precursorList>
        @else:
        @for precursor_element_line in s["precursor_element_string"]:
        @precursor_element_line
        @end
        @end
        @end

        <binaryDataArrayList count="3">
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="mzArray"/>
            <cvParam cvRef="IMS" accession="IMS:1000103" name="external array length" value="@{s["mz_len"]!!s}"/>
            <cvParam cvRef="IMS" accession="IMS:1000104" name="external encoded length" value="@{s["mz_enc_len"]!!s}"/>
            <cvParam cvRef="IMS" accession="IMS:1000102" name="external offset" value="@{s["mz_offset"]!!s}"/>
            <binary/>
          </binaryDataArray>
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="intensityArray"/>
            <cvParam cvRef="IMS" accession="IMS:1000103" name="external array length" value="@{s["int_len"]!!s}"/>
            <cvParam cvRef="IMS" accession="IMS:1000104" name="external encoded length" value="@{s["int_enc_len"]!!s}"/>
            <cvParam cvRef="IMS" accession="IMS:1000102" name="external offset" value="@{s["int_offset"]!!s}"/>
            <binary/>
          </binaryDataArray>
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="mobilityArray"/>
            <cvParam cvRef="IMS" accession="IMS:1000103" name="external array length" value="@{s["mob_len"]!!s}"/>
            <cvParam cvRef="IMS" accession="IMS:1000104" name="external encoded length" value="@{s["mob_enc_len"]!!s}"/>
            <cvParam cvRef="IMS" accession="IMS:1000102" name="external offset" value="@{s["mob_offset"]!!s}"/>
            <binary/>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
      @end
    </spectrumList>
  </run>
</mzML>
"""

class _MaxlenDict(OrderedDict):
    def __init__(self, *args, **kwargs):
        self.maxlen = kwargs.pop('maxlen', None)
        OrderedDict.__init__(self, *args, **kwargs)

    def __setitem__(self, key, value):
        if self.maxlen is not None and len(self) >= self.maxlen:
            self.popitem(0) #pop oldest
        OrderedDict.__setitem__(self, key, value)

class ImzMLWriter(object):
    """
        Create an imzML+ibd file.

        :param output_filename:
            is used to make the base name by removing the extension (if any).
            two files will be made by adding ".ibd" and ".imzML" to the base name
        :param intensity_dtype:
            The numpy data type to use for saving intensity values
        :param mz_dtype:
            The numpy data type to use for saving mz array values
        :param mobility_dtype:
            The numpy data ttype to use for saving mobility array values
        :param mode:

            * "continuous" mode will save the first mz array only
            * "processed" mode save every mz array separately
            * "auto" mode writes only mz arrays that have not already been written
        :param intensity_compression:
            How to compress the intensity data before saving
            must be an instance of :class:`~pyimzml.compression.NoCompression` or :class:`~pyimzml.compression.ZlibCompression`
        :param mz_compression:
            How to compress the mz array data before saving
        :param mobility_compression:
            How to compress the mobility array data before saving
        :param include_mobility:
            bool: True or False
            Whether imzML schema should include trapped ion mobility spectrometry data. Units/metadata based
            on Bruker TIMS data.
        :param polarity:
            str: "positive" or "negative"
            The polarity of the data. If not specified, the polarity will be left blank.
        :param image_x_dimension:
            int: The x dimension of the image in micrometers
        :param image_y_dimension:
            int: The y dimension of the image in micrometers
    """
    def __init__(self, output_filename,
                 mz_dtype=np.float64,
                 intensity_dtype=np.float32,
                 mobility_dtype=np.float64,
                 mode="auto",
                 spec_type="centroid",
                 scan_direction="top_down",
                 line_scan_direction="line_left_right",
                 scan_pattern="flyback",
                 scan_type="horizontal_line",
                 mz_compression=NoCompression(),
                 intensity_compression=NoCompression(),
                 mobility_compression=NoCompression(),
                 polarity=None,
                 include_mobility=False,
                 mobility_info=None,
                 image_x_dimension = None,
                 image_y_dimension = None,
                 xml_element_strings = {}):

        # Whether to include ion mobility data.
        self.include_mobility = include_mobility

        self.ms_levels = []
        self.mz_dtype = mz_dtype
        self.intensity_dtype = intensity_dtype
        self.mobility_dtype = mobility_dtype
        self.mobility_info = mobility_info
        self.mode = mode
        self.spec_type = spec_type
        self.mz_compression = self.compression_string_to_name(mz_compression)
        self.intensity_compression = self.compression_string_to_name(intensity_compression)
        self.mobility_compression = self.compression_string_to_name(mobility_compression)
        self.run_id = os.path.splitext(output_filename)[0]
        self.filename = self.run_id + ".imzML"
        self.ibd_filename = self.run_id + ".ibd"
        self.xml = open(self.filename, 'w')
        self.ibd = open(self.ibd_filename, 'wb+')
        self.sha1 = hashlib.sha1()
        self.uuid = uuid.uuid4()
        
        self.scan_direction = scan_direction
        self.scan_pattern = scan_pattern
        self.scan_type = scan_type
        self.line_scan_direction = line_scan_direction
        self.image_x_dimension = image_x_dimension
        self.image_y_dimension = image_y_dimension
        self.pixel_size_x = None
        self.pixel_size_y = None
        self.xml_element_strings = xml_element_strings

        self._write_ibd(self.uuid.bytes)

        if self.include_mobility == False:
            self.wheezy_engine = Engine(loader=DictLoader({'imzml': IMZML_TEMPLATE}), extensions=[CoreExtension()])
            self.imzml_template = self.wheezy_engine.get_template('imzml')
        elif self.include_mobility == True:
            self.wheezy_engine = Engine(loader=DictLoader({'imzml': IMZML_MOBILITY_TEMPLATE}),
                                        extensions=[CoreExtension()])
            self.imzml_template = self.wheezy_engine.get_template('imzml')
        self.spectra_group = []
        self.spectra = []
        self.first_mz = None
        self.hashes = defaultdict(list)  # mz_hash -> list of mz_data (disk location)
        self.lru_cache = _MaxlenDict(maxlen=10)  # mz_array (as tuple) -> mz_data (disk location)
        self._setPolarity(polarity)

    @staticmethod
    def _np_type_to_name(dtype):
        if dtype.__name__.startswith('float'):
            return "%s-bit float" % dtype.__name__[5:]
        elif dtype.__name__.startswith('int'):
            return "%s-bit integer" % dtype.__name__[3:]
        
    def compression_string_to_name(self, compression_input):
        if compression_input in [NoCompression(), ZlibCompression()]:
            compression_output = compression_input
        elif compression_input is None:
            compression_output = NoCompression()
        elif type(compression_input) == str:
            if compression_input.lower() in ["none", "no compression"]:
                compression_output = NoCompression()
            elif compression_input.lower() in ["zlib", "zlib compression"]:
                compression_output = ZlibCompression()
        else:
            raise ValueError('The input for compression must be "None" or "zlib".')
        return compression_output

    def _setPolarity(self, polarity):
        if polarity:
            if polarity.lower() in ['positive', 'negative']:
                self.polarity = polarity.lower()
            else:
                raise ValueError('value for polarity must be one of "positive", "negative". Received: {}'.format(polarity))
        else:
            self.polarity = ""

    def _write_xml(self):
        spectra = self.spectra
        mz_data_type = self._np_type_to_name(self.mz_dtype)
        int_data_type = self._np_type_to_name(self.intensity_dtype)
        if self.include_mobility == True:
            mob_data_type = self._np_type_to_name(self.mobility_dtype)
            mobility_name, mobility_accession, mobility_unit, mobility_unit_accession = self.mobility_info
        obo_codes = {"32-bit integer": "1000519", 
                     "16-bit float": "1000520",
                     "32-bit float": "1000521",
                     "64-bit integer": "1000522",
                     "64-bit float": "1000523",
                     "continuous": "1000030",
                     "processed": "1000031",
                     "zlib compression": "1000574",
                     "no compression": "1000576",
                     "line_bottom_up": "1000492",
                     "line_left_right": "1000491",
                     "line_right_left": "1000490",
                     "line_top_down": "1000493",
                     "bottom_up": "1000400",
                     "left_right": "1000402",
                     "right_left": "1000403",
                     "top_down": "1000401",
                     "meandering": "1000410",
                     "flyback": "1000413",
                     "random_access": "1000412",
                     "horizontal_line": "1000480",
                     "vertical_line": "1000481"}
        obo_names = {"line_bottom_up": "linescan bottom up",
                     "line_left_right": "linescan left right",
                     "line_right_left": "linescan right left",
                     "line_top_down": "linescan top down",
                     "bottom_up": "bottom up",
                     "left_right": "left right",
                     "right_left": "right left",
                     "top_down": "top down",
                     "meandering": "meandering",
                      "flyback": "flyback",
                     "random_access": "random access",
                     "horizontal_line": "horizontal line scan",
                     "vertical_line": "vertical line scan"}
        
        uuid = ("{%s}" % self.uuid).upper()
        sha1sum = self.sha1.hexdigest().upper()
        run_id = self.run_id
        if self.mode == 'auto':
            mode = "processed" if len(self.lru_cache) > 1 else "continuous"
        else:
            mode = self.mode
        spec_type = self.spec_type
        mz_compression = self.mz_compression.name
        int_compression = self.intensity_compression.name
        if self.include_mobility == True:
            mob_compression = self.mobility_compression.name
        polarity = self.polarity
        scan_direction = self.scan_direction
        scan_pattern = self.scan_pattern
        scan_type = self.scan_type
        line_scan_direction = self.line_scan_direction

        ms_levels = [s["ms_level"] for s in spectra]
        self.ms_levels = []
        for ms_level in ms_levels:
            if ms_level not in self.ms_levels:
                self.ms_levels.append(ms_level)
        ms_levels = self.ms_levels

        image_x_dimension = self.image_x_dimension
        image_y_dimension = self.image_y_dimension
        if image_x_dimension is not None:
            pixel_size_x = self.image_x_dimension/max(s["coords"][0] for s in spectra)
        else:
            pixel_size_x = None
        if image_y_dimension is not None:
            pixel_size_y = self.image_y_dimension/max(s["coords"][1] for s in spectra)
        else:
            pixel_size_y = None
        xml_element_strings = self.xml_element_strings
        self.spectra = spectra
        self.xml.write(self.imzml_template.render(locals()))

    def _write_ibd(self, bytes):
        self.ibd.write(bytes)
        self.sha1.update(bytes)
        return len(bytes)

    def _encode_and_write(self, data, dtype=np.float32, compression=NoCompression()):
        data = np.asarray(data, dtype=dtype)
        offset = self.ibd.tell()
        bytes = data.tobytes()
        bytes = compression.compress(bytes)
        return offset, data.shape[0], self._write_ibd(bytes)

    def _read_mz(self, mz_offset, mz_len, mz_enc_len):
        '''reads a mz array from the currently open ibd file'''
        self.ibd.seek(mz_offset)
        data = self.ibd.read(mz_enc_len)
        self.ibd.seek(0, 2)
        data = self.mz_compression.decompress(data)
        return tuple(np.fromstring(data, dtype=self.mz_dtype))

    def _get_previous_mz(self, mzs):
        '''given an mz array, return the mz_data (disk location)
        if the mz array was not previously written, write to disk first'''
        mzs = tuple(mzs)  # must be hashable
        if mzs in self.lru_cache:
            return self.lru_cache[mzs]

        # mz not recognized ... check hash
        mz_hash = "%s-%s-%s" % (hash(mzs), sum(mzs), len(mzs))
        if mz_hash in self.hashes:
            for mz_data in self.hashes[mz_hash]:
                test_mz = self._read_mz(*mz_data)
                if mzs == test_mz:
                    self.lru_cache[test_mz] = mz_data
                    return mz_data
        # hash not recognized
        # must be a new mz array ... write it, add it to lru_cache and hashes
        mz_data = self._encode_and_write(mzs, self.mz_dtype, self.mz_compression)
        self.hashes[mz_hash].append(mz_data)
        self.lru_cache[mzs] = mz_data
        return mz_data

    def addSpectrum(self, mzs, intensities, coords, mobilities=None, precursor_mz = None, 
                    scan_start_time = None, ms_level = None, filter_string = None, 
                    isolation_window_offset = None, activation = None, mass_window = None,
                    precursor_element_string = None, userParams=[]):
        """
        Add a mass spectrum to the file.

        :param mz:
            mz array
        :param intensities:
            intensity array
        :param mobilities:
            mobility array
        :param coords:

            * 2-tuple of x and y position OR
            * 3-tuple of x, y, and z position

            note some applications want coords to be 1-indexed

        :param precursor_mz:
            mz float
        :param scan_start_time:
            float
        :param ms_level:
            int
        :param filter_string:
            str
        :param isolation_window_offset:
            * float OR
            * 2-tuple of (lower offset, upper offest) as floats
        :param activation:
            Not implemented
        :param mass_window:
            The lower and upper m/z values of the scan window
        """
        # must be rounded now to allow comparisons to later data
        # but don't waste CPU time in continuous mode since the data will not be used anyway
        if self.mode != "continuous" or self.first_mz is None:
            mzs = self.mz_compression.rounding(mzs)
        intensities = self.intensity_compression.rounding(intensities)
        if self.include_mobility == False:
            mobilities = self.mobility_compression.rounding(mobilities)

        if self.mode == "continuous":
            if self.first_mz is None:
                self.first_mz = self._encode_and_write(mzs, self.mz_dtype, self.mz_compression)
            mz_data = self.first_mz
        elif self.mode == "processed":
            mz_data = self._encode_and_write(mzs, self.mz_dtype, self.mz_compression)
        elif self.mode == "auto":
            mz_data = self._get_previous_mz(mzs)
        else:
            raise TypeError("Unknown mode: %s" % self.mode)
        mz_offset, mz_len, mz_enc_len = mz_data

        int_offset, int_len, int_enc_len = self._encode_and_write(intensities, self.intensity_dtype, self.intensity_compression)
        if self.include_mobility == True:
            mob_offset, mob_len, mob_enc_len = self._encode_and_write(mobilities, self.mobility_dtype, self.mobility_compression)
        
        mz_min = np.min(mzs)
        mz_max = np.max(mzs)
        ix_max = np.argmax(intensities)
        mz_base = mzs[ix_max]
        int_base = intensities[ix_max]
        int_tic = np.sum(intensities)
        
        if mass_window is None:
            mass_window = (mz_min, mz_max)

        # Currently ms_level translates to MS1 for level==1 and MSn for any other level.
        # TODO: add support for correct labelling of other scan types such as MRM.
        if ms_level is None:
            if precursor_mz:
                ms_level = 2
            else:
                ms_level = 1
        s = dict(coords=coords, mz_len=mz_len, mz_offset=mz_offset, mz_enc_len=mz_enc_len, 
                 int_len=int_len, int_offset=int_offset, int_enc_len=int_enc_len, 
                 mz_min=mz_min, mz_max=mz_max, mz_base=mz_base, 
                 int_base=int_base, int_tic=int_tic,  activation=activation,
                 scan_start_time=scan_start_time, ms_level=ms_level,
                 filter_string=filter_string, mass_window=mass_window, userParams=userParams)
        if precursor_mz:
            if isolation_window_offset is None:
                isolation_window_lower_offset, isolation_window_upper_offset = None, None
            # if only one offset, use it for both lower and upper, otherwise assume it us lower and upper 
            elif isolation_window_offset is not None:
                try:
                    isolation_window_lower_offset, isolation_window_upper_offset = isolation_window_offset
                except TypeError:
                    isolation_window_lower_offset, isolation_window_upper_offset = isolation_window_offset, isolation_window_offset            
            s.update(precursor_mz=precursor_mz,
                    isolation_window_lower_offset=isolation_window_lower_offset,
                    isolation_window_upper_offset=isolation_window_upper_offset,
                    precursor_element_string=precursor_element_string)
        if self.include_mobility == True:
            s.update(mob_len=mob_len, mob_offset=mob_offset, mob_enc_len=mob_enc_len,
                     mob_min=np.min(mobilities), mob_max=np.max(mobilities))
        self.spectra.append(s)

    def close(self):  # 'close' is a more common use for this
        """
        Writes the XML file and closes all files.
        Will be called automatically if ``with``-pattern is used.
        """
        self.finish()

    def finish(self):
        '''alias of close()'''
        self.ibd.close()
        self._write_xml()
        self.xml.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_t, exc_v, trace):
        if exc_t is None:
            self.finish()
        else:
            self.ibd.close()
            self.xml.close()

# def _main(argv):
#     from pyimzml.ImzMLParser import ImzMLParser
#     inputfile = ''
#     outputfile = ''
#     try:
#         opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
#     except getopt.GetoptError:
#         print('test.py -i <inputfile> -o <outputfile>')
#         sys.exit(2)
#     for opt, arg in opts:
#         if opt == '-h':
#             print('test.py -i <inputfile> -o <outputfile>')
#             sys.exit()
#         elif opt in ("-i", "--ifile"):
#             inputfile = arg
#         elif opt in ("-o", "--ofile"):
#             outputfile = arg
#     if inputfile == '':
#         print('test.py -i <inputfile> -o <outputfile>')
#         raise IOError('input file not specified')
#     if outputfile=='':
#         outputfile=inputfile+'.imzML'
#     imzml = ImzMLParser(inputfile, include_mobility=False)
#     spectra = []
#     with ImzMLWriter(outputfile, mz_dtype=np.float32, intensity_dtype=np.float32, include_mobility=False) as writer:
#         for i, coords in enumerate(imzml.coordinates):
#             mzs, intensities = imzml.getspectrum(i)
#             writer.addSpectrum(mzs, intensities, coords)
#             spectra.append((mzs, intensities, coords))

#     imzml = ImzMLParser(outputfile, include_mobility=False)
#     spectra2 = []
#     for i, coords in enumerate(imzml.coordinates):
#         mzs, intensities = imzml.getspectrum(i)
#         spectra2.append((mzs, intensities, coords))

#     print(spectra[0] == spectra2[0])

#     imzml = ImzMLParser(inputfile, include_mobility=True)
#     spectra = []
#     with ImzMLWriter(outputfile, mz_dtype=np.float32, intensity_dtype=np.float32, include_mobility=True) as writer:
#         for i, coords in enumerate(imzml.coordinates):
#             mzs, intensities, mobilities = imzml.getspectrum(i)
#             writer.addSpectrum(mzs, intensities, coords, mobilities=mobilities)
#             spectra.append((mzs, intensities, mobilities, coords))

#     imzml = ImzMLParser(outputfile, include_mobility=True)
#     spectra2 = []
#     for i, coords in enumerate(imzml.coordinates):
#         mzs, intensities, mobilities = imzml.getspectrum(i)
#         spectra2.append((mzs, intensities, mobilities, coords))

#     print(spectra[0] == spectra2[0])

# if __name__ == '__main__':
#     _main(sys.argv[1:])
