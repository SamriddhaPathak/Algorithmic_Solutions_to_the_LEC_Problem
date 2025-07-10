import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.spatial import SphericalVoronoi
import random

class SphericalVoronoiApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Spherical Voronoi Diagram Generator")
        self.root.geometry("900x700")
        
        # Variables
        self.num_points = tk.IntVar(value=10)
        self.seed_points = []
        self.voronoi = None
        
        self.setup_gui()
        
    def setup_gui(self):
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(1, weight=1)
        
        # Control panel
        control_frame = ttk.LabelFrame(main_frame, text="Controls", padding="10")
        control_frame.grid(row=0, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
        
        # Number of points input
        ttk.Label(control_frame, text="Number of seed points:").grid(row=0, column=0, padx=(0, 10))
        points_spinbox = ttk.Spinbox(control_frame, from_=3, to=100, textvariable=self.num_points, width=10)
        points_spinbox.grid(row=0, column=1, padx=(0, 20))
        
        # Buttons
        ttk.Button(control_frame, text="Generate Random Points", 
                  command=self.generate_random_points).grid(row=0, column=2, padx=(0, 10))
        ttk.Button(control_frame, text="Plot Voronoi Diagram", 
                  command=self.plot_voronoi).grid(row=0, column=3, padx=(0, 10))
        ttk.Button(control_frame, text="Clear", 
                  command=self.clear_plot).grid(row=0, column=4)
        
        # Seed points display
        points_frame = ttk.LabelFrame(main_frame, text="Seed Points (Lat, Lon)", padding="10")
        points_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(0, 10))
        
        # Listbox with scrollbar
        listbox_frame = ttk.Frame(points_frame)
        listbox_frame.pack(fill=tk.BOTH, expand=True)
        
        self.points_listbox = tk.Listbox(listbox_frame, height=15)
        scrollbar = ttk.Scrollbar(listbox_frame, orient=tk.VERTICAL, command=self.points_listbox.yview)
        self.points_listbox.configure(yscrollcommand=scrollbar.set)
        
        self.points_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Add/Remove points buttons
        button_frame = ttk.Frame(points_frame)
        button_frame.pack(fill=tk.X, pady=(10, 0))
        
        ttk.Button(button_frame, text="Add Point", command=self.add_point_dialog).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(button_frame, text="Remove Selected", command=self.remove_point).pack(side=tk.LEFT)
        
        # Plot area
        plot_frame = ttk.LabelFrame(main_frame, text="Voronoi Diagram", padding="10")
        plot_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Create matplotlib figure
        self.fig = plt.Figure(figsize=(8, 6))
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        self.canvas = FigureCanvasTkAgg(self.fig, plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Initialize with empty sphere
        self.draw_empty_sphere()
        
    def generate_random_points(self):
        """Generate random points on the sphere surface"""
        n = self.num_points.get()
        
        # Generate random points on unit sphere using spherical coordinates
        self.seed_points = []
        
        for _ in range(n):
            # Random spherical coordinates
            theta = np.random.uniform(0, 2 * np.pi)  # azimuth
            phi = np.arccos(2 * np.random.uniform(0, 1) - 1)  # polar angle (uniform on sphere)
            
            # Convert to Cartesian coordinates
            x = np.sin(phi) * np.cos(theta)
            y = np.sin(phi) * np.sin(theta)
            z = np.cos(phi)
            
            # Convert to lat/lon for display
            lat = np.degrees(np.pi/2 - phi)
            lon = np.degrees(theta)
            if lon > 180:
                lon -= 360
                
            self.seed_points.append((lat, lon, x, y, z))
        
        self.update_points_listbox()
        
    def add_point_dialog(self):
        """Open dialog to manually add a point"""
        dialog = tk.Toplevel(self.root)
        dialog.title("Add Seed Point")
        dialog.geometry("300x150")
        dialog.transient(self.root)
        dialog.grab_set()
        
        # Center the dialog
        dialog.geometry("+%d+%d" % (self.root.winfo_rootx() + 50, self.root.winfo_rooty() + 50))
        
        # Input fields
        ttk.Label(dialog, text="Latitude (-90 to 90):").pack(pady=5)
        lat_var = tk.DoubleVar(value=0.0)
        lat_entry = ttk.Entry(dialog, textvariable=lat_var)
        lat_entry.pack(pady=5)
        
        ttk.Label(dialog, text="Longitude (-180 to 180):").pack(pady=5)
        lon_var = tk.DoubleVar(value=0.0)
        lon_entry = ttk.Entry(dialog, textvariable=lon_var)
        lon_entry.pack(pady=5)
        
        def add_point():
            try:
                lat = lat_var.get()
                lon = lon_var.get()
                
                if not (-90 <= lat <= 90):
                    raise ValueError("Latitude must be between -90 and 90")
                if not (-180 <= lon <= 180):
                    raise ValueError("Longitude must be between -180 and 180")
                
                # Convert to Cartesian coordinates
                lat_rad = np.radians(lat)
                lon_rad = np.radians(lon)
                
                x = np.cos(lat_rad) * np.cos(lon_rad)
                y = np.cos(lat_rad) * np.sin(lon_rad)
                z = np.sin(lat_rad)
                
                self.seed_points.append((lat, lon, x, y, z))
                self.update_points_listbox()
                dialog.destroy()
                
            except ValueError as e:
                messagebox.showerror("Invalid Input", str(e))
        
        ttk.Button(dialog, text="Add Point", command=add_point).pack(pady=10)
        
    def remove_point(self):
        """Remove selected point from the list"""
        selection = self.points_listbox.curselection()
        if selection:
            index = selection[0]
            del self.seed_points[index]
            self.update_points_listbox()
            
    def update_points_listbox(self):
        """Update the points listbox display"""
        self.points_listbox.delete(0, tk.END)
        for i, (lat, lon, x, y, z) in enumerate(self.seed_points):
            self.points_listbox.insert(tk.END, f"{i+1}: ({lat:.2f}°, {lon:.2f}°)")
            
    def draw_empty_sphere(self):
        """Draw an empty sphere"""
        self.ax.clear()
        
        # Create sphere surface
        u = np.linspace(0, 2 * np.pi, 50)
        v = np.linspace(0, np.pi, 50)
        x_sphere = np.outer(np.cos(u), np.sin(v))
        y_sphere = np.outer(np.sin(u), np.sin(v))
        z_sphere = np.outer(np.ones(np.size(u)), np.cos(v))
        
        # Plot wireframe sphere
        self.ax.plot_wireframe(x_sphere, y_sphere, z_sphere, alpha=0.3, color='lightgray')
        
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        self.ax.set_title('Spherical Surface')
        
        # Set equal aspect ratio
        self.ax.set_box_aspect([1,1,1])
        
        self.canvas.draw()
        
    def plot_voronoi(self):
        """Plot the Voronoi diagram on the sphere"""
        if len(self.seed_points) < 4:
            messagebox.showerror("Error", "Need at least 4 seed points to generate Voronoi diagram")
            return
            
        try:
            # Extract Cartesian coordinates
            points = np.array([[x, y, z] for _, _, x, y, z in self.seed_points])
            
            # Generate spherical Voronoi diagram
            sv = SphericalVoronoi(points, radius=1, center=np.array([0, 0, 0]))
            sv.sort_vertices_of_regions()
            
            self.ax.clear()
            
            # Draw sphere wireframe
            u = np.linspace(0, 2 * np.pi, 30)
            v = np.linspace(0, np.pi, 30)
            x_sphere = np.outer(np.cos(u), np.sin(v))
            y_sphere = np.outer(np.sin(u), np.sin(v))
            z_sphere = np.outer(np.ones(np.size(u)), np.cos(v))
            self.ax.plot_wireframe(x_sphere, y_sphere, z_sphere, alpha=0.1, color='lightgray')
            
            # Plot seed points
            self.ax.scatter(points[:, 0], points[:, 1], points[:, 2], 
                          c='red', s=100, alpha=0.8, label='Seed Points')
            
            # Plot Voronoi vertices
            if hasattr(sv, 'vertices') and len(sv.vertices) > 0:
                self.ax.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2], 
                              c='blue', s=50, alpha=0.6, label='Voronoi Vertices')
            
            # Plot Voronoi edges
            colors = plt.cm.Set3(np.linspace(0, 1, len(sv.regions)))
            
            for region, color in zip(sv.regions, colors):
                if len(region) > 0:
                    # Get vertices of this region
                    region_vertices = sv.vertices[region]
                    
                    # Plot the region boundary
                    for i in range(len(region_vertices)):
                        start = region_vertices[i]
                        end = region_vertices[(i + 1) % len(region_vertices)]
                        self.ax.plot([start[0], end[0]], 
                                   [start[1], end[1]], 
                                   [start[2], end[2]], 
                                   color=color, linewidth=2, alpha=0.8)
            
            self.ax.set_xlabel('X')
            self.ax.set_ylabel('Y')
            self.ax.set_zlabel('Z')
            self.ax.set_title(f'Spherical Voronoi Diagram ({len(self.seed_points)} points)')
            self.ax.legend()
            
            # Set equal aspect ratio
            self.ax.set_box_aspect([1,1,1])
            
            self.canvas.draw()
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate Voronoi diagram: {str(e)}")
            
    def clear_plot(self):
        """Clear the plot and reset"""
        self.seed_points = []
        self.update_points_listbox()
        self.draw_empty_sphere()

def main():
    root = tk.Tk()
    app = SphericalVoronoiApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()