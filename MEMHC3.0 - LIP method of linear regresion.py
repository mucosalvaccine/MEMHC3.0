import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText
import pandas as pd
from mhcflurry import Class1AffinityPredictor
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import threading
import time

# Global dictionary to store data
data_store = {}


def log_message(message):
    log_text.insert(tk.END, message + "\n")
    log_text.see(tk.END)


def generate_peptides(protein_sequence, min_length, max_length):
    peptides = []
    start_positions = []
    for i in range(len(protein_sequence)):
        for j in range(i + min_length, min(i + max_length + 1, len(protein_sequence) + 1)):
            peptides.append(protein_sequence[i:j])
            start_positions.append(i)
    return peptides, start_positions


# Global variables for cumulative data
cumulative_x_data = []
cumulative_y_data = []


def update_chart(x_data, y_data):
    global fig, canvas, ax, cumulative_x_data, cumulative_y_data

    ax.clear()  # Clear the previous chart
    cumulative_x_data.extend(x_data)
    cumulative_y_data.extend(y_data)

    # Plot the data and fill the area under the curve
    ax.plot(cumulative_x_data, cumulative_y_data, color='blue', marker='o', linestyle='-', linewidth=2)
    ax.fill_between(cumulative_x_data, cumulative_y_data, color='blue', alpha=0.7)

    ax.set_xlabel('Percentage of Protein Area Covered')
    ax.set_ylabel('Percentage of HLA Coverage')
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_title("Protein vs HLA Coverage")

    canvas.draw()


def display_protein_matrix(protein_sequence, covered_positions):
    # Create a binary array to represent the coverage of each position
    coverage_array = [1 if i in covered_positions else 0 for i in range(len(protein_sequence))]

    # Clear the previous plot for the matrix
    matrix_ax.clear()

    # Create the matrix heatmap
    matrix_ax.imshow([coverage_array], cmap='Blues', aspect='auto', interpolation='nearest')

    # Set axis labels and ticks
    matrix_ax.set_xticks(range(len(protein_sequence)))
    matrix_ax.set_xticklabels([aa for aa in protein_sequence])
    matrix_ax.set_yticks([])  # No y-axis labels for a single row matrix

    matrix_ax.set_title("Peptide Coverage in Protein Sequence")

    # Redraw the canvas
    matrix_canvas.draw()


def run_script():
    global data_store, cumulative_x_data, cumulative_y_data
    cumulative_x_data.clear()
    cumulative_y_data.clear()

    log_text.delete(1.0, tk.END)
    log_message("Script started.")

    # Clear the treeview
    for row in tree.get_children():
        tree.delete(row)

    protein_sequence = protein_entry.get().upper()
    min_length = int(min_length_entry.get())
    max_length = int(max_length_entry.get())
    hla_file = hla_file_entry.get()
    kd_threshold = int(kd_threshold_entry.get())
    hla_coverage_limit = int(hla_coverage_limit_entry.get())

    if not protein_sequence or not hla_file:
        messagebox.showerror("Input Error", "Please provide all required inputs.")
        return

    # Validate protein sequence
    accepted_aas = set("ACDEFGHIKLMNPQRSTVWY")
    if set(protein_sequence) - accepted_aas:
        messagebox.showerror("Input Error", "Protein sequence contains invalid amino acids.")
        return

    # Generate peptides
    log_message("Generating peptides...")
    peptides, start_positions = generate_peptides(protein_sequence, min_length, max_length)

    # Load predictor and HLA alleles
    log_message("Predicting binding affinities...")
    predictor = Class1AffinityPredictor.load()
    HLA_type_I = pd.read_csv(hla_file)
    HLA_type_I_list = list(HLA_type_I.columns)
    total_hla_alleles = len(HLA_type_I_list)

    # Store all predictions for peptides
    binding_predictions = {}
    for mhc_allele in HLA_type_I_list:
        try:
            predicted_affinities = predictor.predict(peptides, [mhc_allele] * len(peptides))  # Predict for all peptides
            for pep, affinity in zip(peptides, predicted_affinities):
                if pep not in binding_predictions:
                    binding_predictions[pep] = {}
                binding_predictions[pep][mhc_allele] = affinity
        except ValueError:
            pass

    covered_positions = set()  # Set to track non-redundant amino acid positions covered by peptides
    total_covered_alleles = set()  # Track unique HLA alleles covered
    hla_coverage_counts = {hla: 0 for hla in HLA_type_I_list}
    cumulative_protein_coverage = 0

    # Process each peptide
    for pep in peptides:
        hla_hits = [mhc_allele for mhc_allele in HLA_type_I_list if
                    binding_predictions[pep].get(mhc_allele, float('inf')) <= kd_threshold]

        # Filter out HLA alleles that have already been covered more than the limit
        new_hits = [hla for hla in hla_hits if hla_coverage_counts[hla] < hla_coverage_limit]

        # If the peptide has new HLA hits, calculate its coverage
        if new_hits:
            # Calculate absolute HLA coverage
            absolute_coverage = len(new_hits) / total_hla_alleles * 100

            # Track the non-redundant amino acid positions covered by this peptide
            start_pos = start_positions[peptides.index(pep)]
            new_positions = set(range(start_pos, start_pos + len(pep)))

            # Calculate unique new positions
            unique_new_positions = new_positions - covered_positions
            covered_positions.update(unique_new_positions)

            # Update protein coverage
            protein_coverage = (len(unique_new_positions) / len(protein_sequence)) * 100
            cumulative_protein_coverage = (len(covered_positions) / len(protein_sequence)) * 100

            # Update HLA allele coverage
            for hla in new_hits:
                hla_coverage_counts[hla] += 1
            total_covered_alleles.update(new_hits)

            # Calculate cumulative HLA coverage
            cumulative_coverage = len(total_covered_alleles) / total_hla_alleles * 100

            # Update result table and chart
            update_result_table(pep, new_hits, absolute_coverage, cumulative_coverage, cumulative_protein_coverage)
            update_chart([cumulative_protein_coverage], [cumulative_coverage])
            display_protein_matrix(protein_sequence, covered_positions)  # Update the matrix

            # Simulate processing delay to visualize real-time updates
            time.sleep(0.1)

    log_message("Script finished successfully.")
    messagebox.showinfo("Success", "The script has been run successfully.")


def update_result_table(pep, hits, absolute_coverage, cumulative_coverage, protein_coverage):
    # Insert a new row into the treeview (result table)
    tree.insert("", "end", values=(
        pep, ", ".join(hits), f"{absolute_coverage:.2f}", f"{cumulative_coverage:.2f}", f"{protein_coverage:.2f}"))


def browse_file():
    filename = filedialog.askopenfilename()
    hla_file_entry.delete(0, tk.END)
    hla_file_entry.insert(0, filename)


def save_file():
    # Retrieve the data from the treeview
    rows = tree.get_children()
    data = []
    for row in rows:
        data.append(tree.item(row)["values"])

    if not data:
        messagebox.showerror("Error", "No data available to save.")
        return

    # Convert the data to a DataFrame
    df_result = pd.DataFrame(data, columns=['Peptide', 'HLA Hits', 'Absolute Coverage (%)', 'Cumulative Coverage (%)',
                                            'Protein Coverage (%)'])

    # Ask where to save the file
    filetypes = [('CSV Files', '*.csv'), ('Excel Files', '*.xlsx')]
    filepath = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=filetypes)
    if not filepath:
        return

    # Save the data
    if filepath.endswith('.csv'):
        df_result.to_csv(filepath, index=False)
    else:
        df_result.to_excel(filepath, index=False)

    messagebox.showinfo("Save Successful", f"File saved as {filepath}")


def toggle_log():
    log_frame.grid_remove() if log_frame.winfo_ismapped() else log_frame.grid(row=14, column=0, columnspan=2,
                                                                              sticky="nsew")
    toggle_log_button.config(text="Hide Log" if log_frame.winfo_ismapped() else "Show Log")


def start_long_running_task():
    task_thread = threading.Thread(target=run_script)
    task_thread.start()


# Tkinter GUI setup
root = tk.Tk()
root.title("Peptide Generator and HLA Affinity Predictor")
root.geometry("1024x800")

# Input fields
tk.Label(root, text="Protein Sequence:").grid(row=0, column=0, sticky="e")
protein_entry = tk.Entry(root, width=50)
protein_entry.grid(row=0, column=1, padx=10, pady=10)

tk.Label(root, text="Min Peptide Length:").grid(row=1, column=0, sticky="e")
min_length_entry = tk.Entry(root)
min_length_entry.grid(row=1, column=1, padx=10, pady=10)
min_length_entry.insert(0, "8")

tk.Label(root, text="Max Peptide Length:").grid(row=2, column=0, sticky="e")
max_length_entry = tk.Entry(root)
max_length_entry.grid(row=2, column=1, padx=10, pady=10)
max_length_entry.insert(0, "9")

tk.Label(root, text="HLA File:").grid(row=3, column=0, sticky="e")
hla_file_entry = tk.Entry(root, width=50)
hla_file_entry.grid(row=3, column=1, padx=10, pady=10)
tk.Button(root, text="Browse", command=browse_file).grid(row=3, column=2, padx=10, pady=10)

tk.Label(root, text="Kd Threshold (nM):").grid(row=4, column=0, sticky="e")
kd_threshold_entry = tk.Entry(root)
kd_threshold_entry.grid(row=4, column=1, padx=10, pady=10)
kd_threshold_entry.insert(0, "500")

tk.Label(root, text="HLA Coverage Limit:").grid(row=5, column=0, sticky="e")
hla_coverage_limit_entry = tk.Entry(root)
hla_coverage_limit_entry.grid(row=5, column=1, padx=10, pady=10)
hla_coverage_limit_entry.insert(0, "1")

minimize_protein_var = tk.IntVar()
minimize_protein_check = tk.Checkbutton(root, text="Minimize Protein Coverage", variable=minimize_protein_var)
minimize_protein_check.grid(row=6, column=1, sticky="w", padx=10, pady=10)

# Buttons
run_button = tk.Button(root, text="Run LIP method of linear regression", command=start_long_running_task)
run_button.grid(row=7, column=0, sticky="w", padx=10, pady=10)

abort_button = tk.Button(root, text="Abort", command=root.quit)
abort_button.grid(row=7, column=1, sticky="e", padx=10, pady=10)

# Progress labels
tk.Label(root, text="Percentage of HLA Coverage:").grid(row=8, column=0, sticky="e")
hla_coverage_text = tk.StringVar()
hla_coverage_label = tk.Label(root, textvariable=hla_coverage_text)
hla_coverage_label.grid(row=8, column=1, padx=10, pady=10)

# Treeview to display the results
tree = ttk.Treeview(root, height=10)
tree["columns"] = ("Peptide", "HLA Hits", "Absolute Coverage (%)", "Cumulative Coverage (%)", "Protein Coverage (%)")
tree.column("#0", width=0, stretch=tk.NO)
tree.column("Peptide", anchor=tk.W, width=150)
tree.column("HLA Hits", anchor=tk.W, width=300)
tree.column("Absolute Coverage (%)", anchor=tk.CENTER, width=150)
tree.column("Cumulative Coverage (%)", anchor=tk.CENTER, width=150)
tree.column("Protein Coverage (%)", anchor=tk.CENTER, width=150)

tree.heading("Peptide", text="Peptide", anchor=tk.W)
tree.heading("HLA Hits", text="HLA Hits", anchor=tk.W)
tree.heading("Absolute Coverage (%)", text="Absolute Coverage (%)", anchor=tk.CENTER)
tree.heading("Cumulative Coverage (%)", text="Cumulative Coverage (%)", anchor=tk.CENTER)
tree.heading("Protein Coverage (%)", text="Protein Coverage (%)", anchor=tk.CENTER)
tree.grid(row=9, column=0, columnspan=2, padx=10, pady=10, sticky="nsew")

# Save button
save_button = tk.Button(root, text="Download", command=save_file)
save_button.grid(row=10, column=1, padx=10, pady=10)

# Log panel
toggle_log_button = tk.Button(root, text="Show Log", command=toggle_log)
toggle_log_button.grid(row=11, column=1, padx=10, pady=10)

log_frame = tk.Frame(root)
log_text = ScrolledText(log_frame, height=10, state='normal')
log_text.pack(fill=tk.BOTH, expand=True)

root.grid_rowconfigure(9, weight=1)
root.grid_columnconfigure(1, weight=1)

# Live updating chart for HLA coverage vs Protein Area
fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=2, column=3, rowspan=6, padx=10, pady=10, sticky="nsew")

# Add another figure for the protein position matrix
matrix_fig, matrix_ax = plt.subplots(figsize=(10, 1), dpi=100)
matrix_canvas = FigureCanvasTkAgg(matrix_fig, master=root)
matrix_canvas.get_tk_widget().grid(row=8, column=3, rowspan=2, padx=10, pady=10, sticky="nsew")

# Start Tkinter event loop
root.mainloop()
