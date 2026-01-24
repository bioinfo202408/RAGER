# **RAGER table of contents**
1. [Quick start](https://github.com/bioinfo202408/RAGER/blob/main/README.md#quick-start)
2. [Preprocess RNAseq data](https://github.com/bioinfo202408/RAGER/blob/main/Preprocess_RNAseq_data.md)
3. [Preprocess ATACseq data](https://github.com/bioinfo202408/RAGER/blob/main/Preprocess_ATACseq_data.md) 
4. [Joint analysis](https://github.com/bioinfo202408/RAGER/blob/main/Joint_analysis.md)
5. [Custom analysis](https://github.com/bioinfo202408/RAGER/blob/main/Custom_analysis.md)
6. [UI](https://github.com/bioinfo202408/RAGER/blob/main/RAGER_UI.md)

## **RAGER UI README**

This document introduces the **RAGER Streamlit UI** and explains how to launch and use it to run the RAGER Snakemake workflows and download results.

---

## **1. Enter the species-specific directory**

Before starting the UI, please change into the correct working directory based on your species:

```bash
# mouse
cd ~/PROJECT/RAGER/mouse  # Change the directory to PROJECT/RAGER/mouse as the current directory

# human
cd ~/PROJECT/RAGER/human  # Change the directory to PROJECT/RAGER/human as the current directory

# plant
cd ~/PROJECT/RAGER/plant  # Change the directory to PROJECT/RAGER/plant as the current directory
```

---

## **2. Launch the Streamlit UI**

Run the following command to start the UI:

```bash
streamlit run RAGER_UI.py --server.address 127.0.0.1 --server.port 8501
```

### **Access the UI in your browser**

* If your Linux system has a graphical interface, it will typically open the UI automatically in your default browser.
* If you are working on a remote server via SSH (no GUI) and using a laptop/desktop as the client, you need **SSH port forwarding**.

On your **local laptop terminal**, run:

```bash
ssh -L 8501:127.0.0.1:8501 user@server
```

> Replace `user@server` with your real username and server address (e.g., `alice@192.168.1.10`).

Then open the following URL in your **local browser**:

```text
http://127.0.0.1:8501
```

---

## **3. UI modules (run in order)**

The UI provides **four modules**, and it is recommended to run them in the following order:

1. **Preprocess RNAseq**
2. **Preprocess ATACseq**
3. **Joint analysis**
4. **Custom analysis**

### **Prepare input for Custom analysis (required)**

If you plan to run **Custom analysis**, you must create the directory and upload/prepare your custom gene list file on the server in advance.

Create your custom gene list file:

```bash
# Create the custom genes file

mkdir -p ./datasets/Custom_genes_analysis
vim ./datasets/Custom_genes_analysis/custom_genes.txt
```

> The `custom_genes.txt` file should contain one gene per line (recommended).

---

## **4. How to run a module**

Inside each module page (tab), the workflow configuration is controlled by a `config.yaml`.

1. **Edit YAML**: modify parameters in the YAML editor.
2. **Save**: click **Save** to write the updated YAML back to disk.
3. **Run**: click **Run** to launch the Snakemake workflow.
4. **View logs**: the UI prints the running log in real time, and a full log is available in the **Logs** section.

> Tip: You can click **Reload** to discard current edits and reload the original `config.yaml` from disk.
<img width="988" height="834" alt="image" src="https://github.com/user-attachments/assets/61b02531-234b-431c-ae8f-778f537d23b5" />


---

## **5. Results and downloads**

After the workflow finishes, the UI will automatically scan the output directories and display:

* Result tables (e.g., `.csv`, `.txt`)
* Reports and figures in `.pdf` format (previewed when PDF preview is available)

Each result item includes a **Download** button so you can download outputs directly from the UI.

---

## **Notes and Tips**

* **Configuration details**
  For detailed explanations of how to modify the `config.yaml` file, please refer to the corresponding module documentation:
1. [Preprocess RNAseq data](https://github.com/bioinfo202408/RAGER/blob/main/Preprocess_RNAseq_data.md)
2. [Preprocess ATACseq data](https://github.com/bioinfo202408/RAGER/blob/main/Preprocess_ATACseq_data.md) 
3. [Joint analysis](https://github.com/bioinfo202408/RAGER/blob/main/Joint_analysis.md)
4. [Custom analysis](https://github.com/bioinfo202408/RAGER/blob/main/Custom_analysis.md)

Each module has its own Markdown file describing the required parameters and recommended settings.

* **Conda environment**
  Before running any analysis or launching the UI, make sure to activate the RAGER conda environment:

  ```bash
  conda activate rager
  ```

  This environment provides all required dependencies, including Snakemake, Streamlit, and bioinformatics tools used in the workflows.

---






