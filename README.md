# Virtual_Screening_Vina
<p>This pipeline is used for virtual screening using AutoDock Vina using python code.</p>
<p>Contents include how to download the zinc20 database and process it.</p>
<h3>Download zinc20</h3>
<p><code>wget https://files.docking.org/zinc20-ML/ZINC20-ML_smiles.tar.gz</code></p>
<p><code>tar -xvzf ZINC20-ML_smiles.tar.gz
</code></p>
<h3>Processing the zinc20 database</h3>
<ol>
  <li>remove duplicates</li>
  <li>check validity of SMILES strings</li>
  <li>check druglikeness of SMILES strings based on these two following criteria:
      <ul>
        <li>Lipinski’s rule of 5 </li>
        <li>Quantitative Estimate of Druglikeness (QED score) with a cutoff of 0.5</li>
      </ul>
  </li>
</ol>
<p>Run script:
<code>python utils/zinc_clean.py -i &ltinput_zinc_file&gt -o &ltoutput_zinc_file&gt</code></p>
<p><b>Note: Make sure you've already install rdkit package: <code>pip install rdkit</b></code></p>

<h3>Convert smiles to 3D SDF file</h3>
<p>Run code:
<code>python utils/zinc_sdf_prepare.py &ltinput_zinc_file&gt -o &ltoutput_sdf_file&gt</code></p>

