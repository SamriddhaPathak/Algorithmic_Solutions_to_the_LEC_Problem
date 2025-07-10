import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.patches import Circle
from scipy.spatial import Voronoi, voronoi_plot_2d, ConvexHull
import matplotlib.colors as mcolors

class VoronoiLECPlotter:
    def __init__(self, root):
        self.root = root
        self.root.title("Voronoi Diagram & Largest Empty Circle Solver")
        self.root.geometry("1400x900")
        self.root.minsize(1200, 800)
        
        # Configure style
        self.style = ttk.Style()
        self.style.theme_use('clam')
        self.style.configure('TButton', font=('Helvetica', 11), background='#4a86e8', foreground='white')
        self.style.configure('TFrame', background='#f5f5f5')
        self.style.configure('TLabel', font=('Helvetica', 11), background='#f5f5f5')
        
        # Create color scheme for sectors
        self.colors = list(mcolors.TABLEAU_COLORS.values())
        
        # Set up variables
        self.points = []
        self.voronoi_result = None
        self.largest_empty_circle = None
        self.convex_hull = None
        self.show_convex_hull = tk.BooleanVar(value=False)
        
        # Define original axis limits
        self.original_xlim = [0, 10]
        self.original_ylim = [0, 10]
        
        # Create UI elements
        self.create_ui()
        
        # Add help information to status bar on startup
        self.status_var.set("Ready. Click to add sectors. Use toolbar to zoom/pan. Home button to reset view.")
        
    def create_ui(self):
        # Main frame
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Left panel for controls
        control_frame = ttk.Frame(main_frame, width=300)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10), pady=10)
        
        # Title
        title_label = ttk.Label(control_frame, text="Voronoi & LEC Solver", font=('Helvetica', 14, 'bold'))
        title_label.pack(pady=(0, 15))
        
        # Instructions
        instruction_label = ttk.Label(control_frame, text="Click on the right panel to add sectors or use manual entry below.\nThen click 'Generate' to create the diagram and find the largest empty circle.", wraplength=280)
        instruction_label.pack(pady=(0, 15))
        
        # Manual entry section
        manual_frame = ttk.LabelFrame(control_frame, text="Manual Sector Entry")
        manual_frame.pack(fill=tk.X, pady=(0, 15))
        
        # X coordinate entry
        x_frame = ttk.Frame(manual_frame)
        x_frame.pack(fill=tk.X, padx=5, pady=5)
        ttk.Label(x_frame, text="X:").pack(side=tk.LEFT)
        self.x_entry = ttk.Entry(x_frame, width=10)
        self.x_entry.pack(side=tk.LEFT, padx=(5, 0))
        
        # Y coordinate entry
        y_frame = ttk.Frame(manual_frame)
        y_frame.pack(fill=tk.X, padx=5, pady=5)
        ttk.Label(y_frame, text="Y:").pack(side=tk.LEFT)
        self.y_entry = ttk.Entry(y_frame, width=10)
        self.y_entry.pack(side=tk.LEFT, padx=(5, 0))
        
        # Add manual point button
        add_manual_button = ttk.Button(manual_frame, text="Add Sector", command=self.add_manual_point)
        add_manual_button.pack(fill=tk.X, padx=5, pady=5)
        
        # Bind Enter key to add manual point
        self.x_entry.bind('<Return>', lambda e: self.add_manual_point())
        self.y_entry.bind('<Return>', lambda e: self.add_manual_point())
        
        # Options section
        options_frame = ttk.LabelFrame(control_frame, text="Display Options")
        options_frame.pack(fill=tk.X, pady=(0, 15))
        
        # Convex hull checkbox
        self.convex_hull_checkbox = ttk.Checkbutton(options_frame, text="Show Convex Hull", 
                                                   variable=self.show_convex_hull,
                                                   command=self.toggle_convex_hull)
        self.convex_hull_checkbox.pack(anchor=tk.W, padx=5, pady=5)
        
        # Buttons frame
        buttons_frame = ttk.Frame(control_frame)
        buttons_frame.pack(fill=tk.X, pady=10)
        
        # Generate button
        self.generate_button = ttk.Button(buttons_frame, text="Generate Voronoi & LEC", command=self.generate_voronoi_and_lec)
        self.generate_button.pack(fill=tk.X, pady=5)
        
        # Clear button
        self.clear_button = ttk.Button(buttons_frame, text="Clear All", command=self.clear_all)
        self.clear_button.pack(fill=tk.X, pady=5)
        
        # Points list label
        points_label = ttk.Label(control_frame, text="Added Sectors:", font=('Helvetica', 11, 'bold'))
        points_label.pack(pady=(20, 5), anchor=tk.W)
        
        # Points list frame
        points_frame = ttk.Frame(control_frame)
        points_frame.pack(fill=tk.BOTH, expand=True)
        
        # Scrollable points list
        self.points_listbox = tk.Listbox(points_frame, font=('Helvetica', 10), bg='white', 
                                         selectbackground='#4a86e8', selectforeground='white')
        self.points_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # Scrollbar for points list
        points_scrollbar = ttk.Scrollbar(points_frame, orient="vertical", command=self.points_listbox.yview)
        points_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.points_listbox.config(yscrollcommand=points_scrollbar.set)
        
        # Remove selected point button
        self.remove_button = ttk.Button(control_frame, text="Remove Selected", command=self.remove_selected)
        self.remove_button.pack(fill=tk.X, pady=5)
        
        # Right side container
        right_container = ttk.Frame(main_frame)
        right_container.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        # Plot panel
        plot_frame = ttk.Frame(right_container)
        plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        # Create matplotlib figure
        self.fig = plt.Figure(figsize=(6, 5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlim(self.original_xlim)
        self.ax.set_ylim(self.original_ylim)
        self.ax.set_title('Click to Add Sectors')
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.grid(True, linestyle='--', alpha=0.7)
        
        # Embed matplotlib figure in tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)
        
        # Add navigation toolbar for zoom and pan functionality
        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        self.toolbar.update()
        
        # Bind click event
        self.canvas.mpl_connect('button_press_event', self.on_plot_click)
        
        # Information panel
        info_frame = ttk.LabelFrame(right_container, text="Diagram Information")
        info_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=10)
        
        # Create tabs for different information sections
        self.info_tabs = ttk.Notebook(info_frame)
        self.info_tabs.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Sites tab
        sites_frame = ttk.Frame(self.info_tabs)
        self.info_tabs.add(sites_frame, text="Sites")
        
        self.sites_text = scrolledtext.ScrolledText(sites_frame, height=6, wrap=tk.WORD)
        self.sites_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Vertices tab
        vertices_frame = ttk.Frame(self.info_tabs)
        self.info_tabs.add(vertices_frame, text="Vertices")
        
        self.vertices_text = scrolledtext.ScrolledText(vertices_frame, height=6, wrap=tk.WORD)
        self.vertices_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # LEC tab
        lec_frame = ttk.Frame(self.info_tabs)
        self.info_tabs.add(lec_frame, text="Largest Empty Circle")
        
        self.lec_text = scrolledtext.ScrolledText(lec_frame, height=6, wrap=tk.WORD)
        self.lec_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Convex Hull tab
        hull_frame = ttk.Frame(self.info_tabs)
        self.info_tabs.add(hull_frame, text="Convex Hull")
        
        self.hull_text = scrolledtext.ScrolledText(hull_frame, height=6, wrap=tk.WORD)
        self.hull_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Status bar
        self.status_var = tk.StringVar()
        self.status_var.set("Ready. Click on the plot to add sectors.")
        status_bar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)
    
    def add_manual_point(self):
        try:
            x = float(self.x_entry.get())
            y = float(self.y_entry.get())
            
            # Add point to the list
            self.points.append((x, y))
            point_idx = len(self.points)
            
            # Update the list display
            self.points_listbox.insert(tk.END, f"Sector {point_idx}: ({x:.2f}, {y:.2f})")
            
            # Update plot
            color = self.colors[(point_idx - 1) % len(self.colors)]
            self.ax.plot(x, y, 'o', color=color, markersize=10)
            self.ax.text(x+0.1, y+0.1, f"{point_idx}", fontsize=9, weight='bold')
            self.canvas.draw()
            
            # Clear entry fields
            self.x_entry.delete(0, tk.END)
            self.y_entry.delete(0, tk.END)
            self.x_entry.focus()
            
            # Update status
            self.status_var.set(f"Added Sector {point_idx} at ({x:.2f}, {y:.2f})")
            
            # Clear any previous information
            self.sites_text.delete(1.0, tk.END)
            self.vertices_text.delete(1.0, tk.END)
            self.lec_text.delete(1.0, tk.END)
            self.hull_text.delete(1.0, tk.END)
            
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter valid numeric coordinates.")
            self.x_entry.focus()
    
    def on_plot_click(self, event):
        if event.inaxes == self.ax:
            # Check if toolbar is in pan/zoom mode
            if self.toolbar.mode != '':
                return  # Skip point addition if in navigate mode
                
            x, y = event.xdata, event.ydata
            if x is not None and y is not None:
                # Add point to the list
                self.points.append((x, y))
                point_idx = len(self.points)
                
                # Update the list display
                self.points_listbox.insert(tk.END, f"Sector {point_idx}: ({x:.2f}, {y:.2f})")
                
                # Update plot
                color = self.colors[(point_idx - 1) % len(self.colors)]
                self.ax.plot(x, y, 'o', color=color, markersize=10)
                self.ax.text(x+0.1, y+0.1, f"{point_idx}", fontsize=9, weight='bold')
                self.canvas.draw()
                
                # Update status
                self.status_var.set(f"Added Sector {point_idx} at ({x:.2f}, {y:.2f})")
                
                # Clear any previous information
                self.sites_text.delete(1.0, tk.END)
                self.vertices_text.delete(1.0, tk.END)
                self.lec_text.delete(1.0, tk.END)
                self.hull_text.delete(1.0, tk.END)
    
    def compute_convex_hull(self):
        if len(self.points) < 3:
            self.convex_hull = None
            return
        
        points = np.array(self.points)
        try:
            hull = ConvexHull(points)
            self.convex_hull = hull
        except Exception as e:
            print(f"Error computing convex hull: {e}")
            self.convex_hull = None
    
    def draw_convex_hull(self):
        if self.convex_hull is not None and self.show_convex_hull.get():
            points = np.array(self.points)
            hull_points = points[self.convex_hull.vertices]
            
            # Close the hull by adding the first point at the end
            hull_points = np.vstack([hull_points, hull_points[0]])
            
            # Draw the convex hull
            self.ax.plot(hull_points[:, 0], hull_points[:, 1], 'g-', linewidth=2, alpha=0.7, label='Convex Hull')
            self.ax.fill(hull_points[:, 0], hull_points[:, 1], 'green', alpha=0.1)
    
    def toggle_convex_hull(self):
        if len(self.points) >= 3:
            # Redraw the plot with or without convex hull
            self.redraw_plot()
    
    def redraw_plot(self):
        # Get current axes limits before updating
        x_lim = self.ax.get_xlim()
        y_lim = self.ax.get_ylim()
        
        # Clear the plot
        self.ax.clear()
        
        # Reset to the previous view limits
        self.ax.set_xlim(x_lim)
        self.ax.set_ylim(y_lim)
        
        # Set title based on current state
        if self.voronoi_result is not None:
            self.ax.set_title('Voronoi Diagram & Largest Empty Circle')
        else:
            self.ax.set_title('Click to Add Sectors')
            
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.grid(True, linestyle='--', alpha=0.7)
        
        # Draw convex hull first (so it appears behind other elements)
        if len(self.points) >= 3:
            self.compute_convex_hull()
            self.draw_convex_hull()
        
        # Plot the Voronoi diagram if it exists
        if self.voronoi_result is not None:
            voronoi_plot_2d(self.voronoi_result, ax=self.ax, show_vertices=True, 
                         line_colors='black', line_width=2, 
                         line_alpha=0.6, point_size=10)
            
            # Draw largest empty circle if it exists
            if self.largest_empty_circle is not None and self.largest_empty_circle['center'] is not None:
                center = self.largest_empty_circle['center']
                radius = self.largest_empty_circle['radius']
                circle = Circle(center, radius, fill=False, 
                               edgecolor='red', linestyle='--', linewidth=2, alpha=0.8)
                self.ax.add_patch(circle)
                self.ax.plot(center[0], center[1], 'rx', markersize=8)
        
        # Replot the original points with colors and labels
        for i, (x, y) in enumerate(self.points):
            color = self.colors[i % len(self.colors)]
            self.ax.plot(x, y, 'o', color=color, markersize=10)
            self.ax.text(x+0.1, y+0.1, f"{i+1}", fontsize=9, weight='bold')
        
        # Update the canvas
        self.canvas.draw()
    
    def generate_voronoi_and_lec(self):
        if len(self.points) < 3:
            messagebox.showwarning("Not Enough Points", "Please add at least 3 sectors to generate a Voronoi diagram.")
            return
        
        # Convert points list to numpy array
        points = np.array(self.points)
        
        # Generate Voronoi diagram
        vor = Voronoi(points)
        self.voronoi_result = vor
        
        # Find largest empty circle
        self.find_largest_empty_circle(vor)
        
        # Compute convex hull
        self.compute_convex_hull()
        
        # Redraw everything
        self.redraw_plot()
        
        # Update information panels
        self.update_information_panels()
        
        self.status_var.set("Voronoi diagram and Largest Empty Circle generated successfully!")
    
    def find_largest_empty_circle(self, vor):
        # For LEC, we need to consider:
        # 1. All Voronoi vertices (circle centers with 3 sites on boundary)
        # 2. The farthest Voronoi vertex from each input site
        
        points = np.array(self.points)
        vertices = vor.vertices
        
        # Calculate distance from each vertex to the nearest site
        best_radius = 0
        best_center = None
        
        # Check if vertices exist (could be none in degenerate cases)
        if len(vertices) > 0:
            for vertex in vertices:
                # Calculate distances to all sites
                distances = np.sqrt(np.sum((points - vertex)**2, axis=1))
                min_distance = np.min(distances)
                
                # Update best circle if this one is larger
                if min_distance > best_radius:
                    best_radius = min_distance
                    best_center = vertex
        
        # Store the result
        self.largest_empty_circle = {
            'center': best_center,
            'radius': best_radius
        }
    
    def update_information_panels(self):
        # Update Sites information
        self.sites_text.delete(1.0, tk.END)
        self.sites_text.insert(tk.END, "Site Coordinates:\n")
        for i, point in enumerate(self.points):
            self.sites_text.insert(tk.END, f"Site {i+1}: ({point[0]:.4f}, {point[1]:.4f})\n")
        
        # Update Vertices information
        if self.voronoi_result is not None:
            self.vertices_text.delete(1.0, tk.END)
            self.vertices_text.insert(tk.END, "Voronoi Vertices Coordinates:\n")
            for i, vertex in enumerate(self.voronoi_result.vertices):
                self.vertices_text.insert(tk.END, f"Vertex {i+1}: ({vertex[0]:.4f}, {vertex[1]:.4f})\n")
        
        # Update LEC information
        if self.largest_empty_circle is not None:
            self.lec_text.delete(1.0, tk.END)
            center = self.largest_empty_circle['center']
            radius = self.largest_empty_circle['radius']
            self.lec_text.insert(tk.END, "Largest Empty Circle:\n")
            if center is not None:
                self.lec_text.insert(tk.END, f"Center: ({center[0]:.4f}, {center[1]:.4f})\n")
                self.lec_text.insert(tk.END, f"Radius: {radius:.4f}\n")
                
                # Find nearest sites (those that define the circle)
                points = np.array(self.points)
                distances = np.sqrt(np.sum((points - center)**2, axis=1))
                nearest_indices = np.argsort(distances)[:3]
                
                self.lec_text.insert(tk.END, "\nNearest Sites (likely on circle boundary):\n")
                for idx in nearest_indices:
                    distance = distances[idx]
                    if abs(distance - radius) < 0.001:  # Consider it's on boundary
                        self.lec_text.insert(tk.END, f"Site {idx+1}: ({points[idx][0]:.4f}, {points[idx][1]:.4f}) - On boundary\n")
                    else:
                        self.lec_text.insert(tk.END, f"Site {idx+1}: ({points[idx][0]:.4f}, {points[idx][1]:.4f}) - Distance: {distance:.4f}\n")
            else:
                self.lec_text.insert(tk.END, "No valid largest empty circle found.\n")
        
        # Update Convex Hull information
        if self.convex_hull is not None:
            self.hull_text.delete(1.0, tk.END)
            self.hull_text.insert(tk.END, "Convex Hull Information:\n")
            self.hull_text.insert(tk.END, f"Number of vertices: {len(self.convex_hull.vertices)}\n")
            self.hull_text.insert(tk.END, f"Area: {self.convex_hull.volume:.4f}\n\n")
            
            self.hull_text.insert(tk.END, "Hull Vertices (in order):\n")
            points = np.array(self.points)
            for i, vertex_idx in enumerate(self.convex_hull.vertices):
                point = points[vertex_idx]
                self.hull_text.insert(tk.END, f"Vertex {i+1}: Site {vertex_idx+1} at ({point[0]:.4f}, {point[1]:.4f})\n")
        else:
            self.hull_text.delete(1.0, tk.END)
            if len(self.points) < 3:
                self.hull_text.insert(tk.END, "Need at least 3 points to compute convex hull.\n")
            else:
                self.hull_text.insert(tk.END, "Convex hull computation failed.\n")
    
    def clear_all(self):
        # Clear the points list
        self.points = []
        self.points_listbox.delete(0, tk.END)
        
        # Clear entry fields
        self.x_entry.delete(0, tk.END)
        self.y_entry.delete(0, tk.END)
        
        # Clear the plot and reset to original view limits
        self.ax.clear()
        self.ax.set_xlim(self.original_xlim)
        self.ax.set_ylim(self.original_ylim)
        self.ax.set_title('Click to Add Sectors')
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.grid(True, linestyle='--', alpha=0.7)
        self.canvas.draw()
        
        # Clear information panels
        self.sites_text.delete(1.0, tk.END)
        self.vertices_text.delete(1.0, tk.END)
        self.lec_text.delete(1.0, tk.END)
        self.hull_text.delete(1.0, tk.END)
        
        # Reset variables
        self.voronoi_result = None
        self.largest_empty_circle = None
        self.convex_hull = None
        
        # Update status
        self.status_var.set("All sectors cleared. Click on the plot to add new sectors.")
    
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
        # Keep current view if there are still points, otherwise reset to original
        if len(self.points) > 0:
            # Get current axes limits before updating
            x_lim = self.ax.get_xlim()
            y_lim = self.ax.get_ylim()
            self.ax.set_xlim(x_lim)
            self.ax.set_ylim(y_lim)
        else:
            # Reset to original view limits if all points are removed
            self.ax.set_xlim(self.original_xlim)
            self.ax.set_ylim(self.original_ylim)
            
        self.ax.set_title('Click to Add Sectors')
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.grid(True, linestyle='--', alpha=0.7)
        
        # Replot the points
        for i, (x, y) in enumerate(self.points):
            color = self.colors[i % len(self.colors)]
            self.ax.plot(x, y, 'o', color=color, markersize=10)
            self.ax.text(x+0.1, y+0.1, f"{i+1}", fontsize=9, weight='bold')
        
        self.canvas.draw()
        
        # Clear information panels
        self.sites_text.delete(1.0, tk.END)
        self.vertices_text.delete(1.0, tk.END)
        self.lec_text.delete(1.0, tk.END)
        self.hull_text.delete(1.0, tk.END)
        
        # Reset variables
        self.voronoi_result = None
        self.largest_empty_circle = None
        self.convex_hull = None
        
        self.status_var.set(f"Removed sector. {len(self.points)} sectors remaining.")

if __name__ == "__main__":
    root = tk.Tk()
    app = VoronoiLECPlotter(root)
    root.mainloop()