import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.colors as mcolors


class VoronoiPlotter:
    def __init__(self, root):
        self.root = root
        self.root.title("Voronoi Diagram Plotter")
        self.root.geometry("1000x700")
        self.root.minsize(900, 600)

        # Configure style
        self.style = ttk.Style()
        self.style.theme_use("clam")
        self.style.configure(
            "TButton", font=("Helvetica", 11), background="#4a86e8", foreground="white"
        )
        self.style.configure("TFrame", background="#f5f5f5")
        self.style.configure("TLabel", font=("Helvetica", 11), background="#f5f5f5")

        # Create color scheme for sectors
        self.colors = list(mcolors.TABLEAU_COLORS.values())

        # Set up variables
        self.points = []

        # Create UI elements
        self.create_ui()

    def create_ui(self):
        # Main frame
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Left panel for controls
        control_frame = ttk.Frame(main_frame, width=250)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10), pady=10)

        # Title
        title_label = ttk.Label(
            control_frame,
            text="Voronoi Diagram Plotter",
            font=("Helvetica", 14, "bold"),
        )
        title_label.pack(pady=(0, 20))

        # Instructions
        instruction_label = ttk.Label(
            control_frame,
            text="Click on the right panel to add sectors.\nThen click 'Generate Voronoi' to create the diagram.",
            wraplength=230,
        )
        instruction_label.pack(pady=(0, 20))

        # Buttons frame
        buttons_frame = ttk.Frame(control_frame)
        buttons_frame.pack(fill=tk.X, pady=10)

        # Generate button
        self.generate_button = ttk.Button(
            buttons_frame, text="Generate Voronoi", command=self.generate_voronoi
        )
        self.generate_button.pack(fill=tk.X, pady=5)

        # Clear button
        self.clear_button = ttk.Button(
            buttons_frame, text="Clear All", command=self.clear_all
        )
        self.clear_button.pack(fill=tk.X, pady=5)

        # Points list label
        points_label = ttk.Label(
            control_frame, text="Added Sectors:", font=("Helvetica", 11, "bold")
        )
        points_label.pack(pady=(20, 5), anchor=tk.W)

        # Points list frame
        points_frame = ttk.Frame(control_frame)
        points_frame.pack(fill=tk.BOTH, expand=True)

        # Scrollable points list
        self.points_listbox = tk.Listbox(
            points_frame,
            font=("Helvetica", 10),
            bg="white",
            selectbackground="#4a86e8",
            selectforeground="white",
        )
        self.points_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Scrollbar for points list
        points_scrollbar = ttk.Scrollbar(
            points_frame, orient="vertical", command=self.points_listbox.yview
        )
        points_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.points_listbox.config(yscrollcommand=points_scrollbar.set)

        # Remove selected point button
        self.remove_button = ttk.Button(
            control_frame, text="Remove Selected", command=self.remove_selected
        )
        self.remove_button.pack(fill=tk.X, pady=5)

        # Right panel for plot
        plot_frame = ttk.Frame(main_frame)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Create matplotlib figure
        self.fig = plt.Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlim([0, 10])
        self.ax.set_ylim([0, 10])
        self.ax.set_title("Click to Add Sectors")
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.grid(True, linestyle="--", alpha=0.7)

        # Embed matplotlib figure in tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)

        # Bind click event
        self.canvas.mpl_connect("button_press_event", self.on_plot_click)

        # Status bar
        self.status_var = tk.StringVar()
        self.status_var.set("Ready. Click on the plot to add sectors.")
        status_bar = ttk.Label(
            self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W
        )
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)

    def on_plot_click(self, event):
        if event.inaxes == self.ax:
            x, y = event.xdata, event.ydata
            if x is not None and y is not None:
                # Add point to the list
                self.points.append((x, y))
                point_idx = len(self.points)

                # Update the list display
                self.points_listbox.insert(
                    tk.END, f"Sector {point_idx}: ({x:.2f}, {y:.2f})"
                )

                # Update plot
                color = self.colors[(point_idx - 1) % len(self.colors)]
                self.ax.plot(x, y, "o", color=color, markersize=10)
                self.ax.text(
                    x + 0.1, y + 0.1, f"{point_idx}", fontsize=9, weight="bold"
                )
                self.canvas.draw()

                # Update status
                self.status_var.set(f"Added Sector {point_idx} at ({x:.2f}, {y:.2f})")

    def generate_voronoi(self):
        if len(self.points) < 3:
            messagebox.showwarning(
                "Not Enough Points",
                "Please add at least 3 sectors to generate a Voronoi diagram.",
            )
            return

        # Convert points list to numpy array
        points = np.array(self.points)

        # Generate Voronoi diagram
        vor = Voronoi(points)

        # Clear the plot but keep the points
        self.ax.clear()
        self.ax.set_xlim([0, 10])
        self.ax.set_ylim([0, 10])
        self.ax.set_title("Voronoi Diagram")
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.grid(True, linestyle="--", alpha=0.7)

        # Plot the Voronoi diagram
        voronoi_plot_2d(
            vor,
            ax=self.ax,
            show_vertices=False,
            line_colors="black",
            line_width=2,
            line_alpha=0.6,
            point_size=10,
        )

        # Replot the original points with colors and labels
        for i, (x, y) in enumerate(self.points):
            color = self.colors[i % len(self.colors)]
            self.ax.plot(x, y, "o", color=color, markersize=10)
            self.ax.text(x + 0.1, y + 0.1, f"{i+1}", fontsize=9, weight="bold")

        # Update the canvas
        self.canvas.draw()
        self.status_var.set("Voronoi diagram generated successfully!")

    def clear_all(self):
        # Clear the points list
        self.points = []
        self.points_listbox.delete(0, tk.END)

        # Clear the plot
        self.ax.clear()
        self.ax.set_xlim([0, 10])
        self.ax.set_ylim([0, 10])
        self.ax.set_title("Click to Add Sectors")
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.grid(True, linestyle="--", alpha=0.7)
        self.canvas.draw()

        # Update status
        self.status_var.set(
            "All sectors cleared. Click on the plot to add new sectors."
        )

    def remove_selected(self):
        selected = self.points_listbox.curselection()
        if not selected:
            return

        # Get the index of the selected point
        idx = selected[0]

        # Remove from points list
        self.points.pop(idx)

        # Remove from listbox
        self.points_listbox.delete(idx)

        # Update the indices in the listbox
        for i in range(idx, self.points_listbox.size()):
            old_text = self.points_listbox.get(i)
            new_text = f"Sector {i+1}:" + old_text.split(":", 1)[1]
            self.points_listbox.delete(i)
            self.points_listbox.insert(i, new_text)

        # Redraw the plot
        self.ax.clear()
        self.ax.set_xlim([0, 10])
        self.ax.set_ylim([0, 10])
        self.ax.set_title("Click to Add Sectors")
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.grid(True, linestyle="--", alpha=0.7)

        # Replot the points
        for i, (x, y) in enumerate(self.points):
            color = self.colors[i % len(self.colors)]
            self.ax.plot(x, y, "o", color=color, markersize=10)
            self.ax.text(x + 0.1, y + 0.1, f"{i+1}", fontsize=9, weight="bold")

        self.canvas.draw()
        self.status_var.set(f"Removed sector. {len(self.points)} sectors remaining.")


if __name__ == "__main__":
    root = tk.Tk()
    app = VoronoiPlotter(root)
    root.mainloop()
