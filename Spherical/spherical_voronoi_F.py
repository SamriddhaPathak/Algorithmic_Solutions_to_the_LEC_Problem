import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.spatial import SphericalVoronoi
from scipy.optimize import minimize
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
        self.largest_empty_circle = None
        self.empty_circle_radius = 0.0

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
        control_frame.grid(
            row=0, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10)
        )

        # Number of points input
        ttk.Label(control_frame, text="Number of seed points:").grid(
            row=0, column=0, padx=(0, 10)
        )
        points_spinbox = ttk.Spinbox(
            control_frame, from_=3, to=100, textvariable=self.num_points, width=10
        )
        points_spinbox.grid(row=0, column=1, padx=(0, 20))

        # Buttons
        ttk.Button(
            control_frame,
            text="Generate Random Points",
            command=self.generate_random_points,
        ).grid(row=0, column=2, padx=(0, 10))
        ttk.Button(
            control_frame, text="Plot Voronoi Diagram", command=self.plot_voronoi
        ).grid(row=0, column=3, padx=(0, 10))
        ttk.Button(control_frame, text="Clear", command=self.clear_plot).grid(
            row=0, column=4
        )

        # Empty circle info display
        info_frame = ttk.Frame(control_frame)
        info_frame.grid(
            row=1, column=0, columnspan=5, pady=(10, 0), sticky=(tk.W, tk.E)
        )

        self.empty_circle_label = ttk.Label(
            info_frame,
            text="Largest Empty Circle Radius: Not calculated",
            font=("Arial", 10, "bold"),
        )
        self.empty_circle_label.pack(side=tk.LEFT)

        # Seed points display
        points_frame = ttk.LabelFrame(
            main_frame, text="Seed Points (Lat, Lon)", padding="10"
        )
        points_frame.grid(
            row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(0, 10)
        )

        # Listbox with scrollbar
        listbox_frame = ttk.Frame(points_frame)
        listbox_frame.pack(fill=tk.BOTH, expand=True)

        self.points_listbox = tk.Listbox(listbox_frame, height=15)
        scrollbar = ttk.Scrollbar(
            listbox_frame, orient=tk.VERTICAL, command=self.points_listbox.yview
        )
        self.points_listbox.configure(yscrollcommand=scrollbar.set)

        self.points_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Add/Remove points buttons
        button_frame = ttk.Frame(points_frame)
        button_frame.pack(fill=tk.X, pady=(10, 0))

        ttk.Button(button_frame, text="Add Point", command=self.add_point_dialog).pack(
            side=tk.LEFT, padx=(0, 5)
        )
        ttk.Button(
            button_frame, text="Remove Selected", command=self.remove_point
        ).pack(side=tk.LEFT)

        # Plot area
        plot_frame = ttk.LabelFrame(main_frame, text="Voronoi Diagram", padding="10")
        plot_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Create matplotlib figure
        self.fig = plt.Figure(figsize=(8, 6))
        self.ax = self.fig.add_subplot(111, projection="3d")

        self.canvas = FigureCanvasTkAgg(self.fig, plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Initialize with empty sphere
        self.draw_empty_sphere()

    def generate_random_points(self):
        """Generate random points on the sphere surface"""
        n = self.num_points.get()

        # Clear existing points first
        self.seed_points = []

        # Generate random points on unit sphere using spherical coordinates
        np.random.seed()  # Ensure we get different random points each time

        for _ in range(n):
            # Random spherical coordinates
            theta = np.random.uniform(0, 2 * np.pi)  # azimuth
            phi = np.arccos(
                2 * np.random.uniform(0, 1) - 1
            )  # polar angle (uniform on sphere)

            # Convert to Cartesian coordinates
            x = np.sin(phi) * np.cos(theta)
            y = np.sin(phi) * np.sin(theta)
            z = np.cos(phi)

            # Convert to lat/lon for display
            lat = np.degrees(np.pi / 2 - phi)
            lon = np.degrees(theta)
            if lon > 180:
                lon -= 360

            self.seed_points.append((lat, lon, x, y, z))

        # Reset empty circle info
        self.largest_empty_circle = None
        self.empty_circle_radius = 0.0
        self.empty_circle_label.config(
            text="Largest Empty Circle Radius: Not calculated"
        )

        self.update_points_listbox()

        # Display the generated points immediately
        self.display_points_only()

    def add_point_dialog(self):
        """Open dialog to manually add a point"""
        dialog = tk.Toplevel(self.root)
        dialog.title("Add Seed Point")
        dialog.geometry("300x150")
        dialog.transient(self.root)
        dialog.grab_set()

        # Center the dialog
        dialog.geometry(
            "+%d+%d" % (self.root.winfo_rootx() + 50, self.root.winfo_rooty() + 50)
        )

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
                # Reset empty circle info when points change
                self.largest_empty_circle = None
                self.empty_circle_radius = 0.0
                self.empty_circle_label.config(
                    text="Largest Empty Circle Radius: Not calculated"
                )
                self.update_points_listbox()
                # Update display to show the new point
                self.display_points_only()
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
            # Reset empty circle info when points change
            self.largest_empty_circle = None
            self.empty_circle_radius = 0.0
            self.empty_circle_label.config(
                text="Largest Empty Circle Radius: Not calculated"
            )
            self.update_points_listbox()
            # Update display
            self.display_points_only()

    def spherical_distance(self, p1, p2):
        """Calculate the great circle distance between two points on a unit sphere"""
        # Ensure points are normalized (on unit sphere)
        p1 = p1 / np.linalg.norm(p1)
        p2 = p2 / np.linalg.norm(p2)

        # Calculate dot product and clamp to avoid numerical errors
        dot_product = np.clip(np.dot(p1, p2), -1.0, 1.0)

        # Return the angular distance (in radians)
        return np.arccos(dot_product)

    def find_largest_empty_circle(self, points):
        """Find the largest empty circle on the sphere"""

        def objective(sphere_coords):
            # Convert spherical coordinates to Cartesian
            theta, phi = sphere_coords
            x = np.sin(phi) * np.cos(theta)
            y = np.sin(phi) * np.sin(theta)
            z = np.cos(phi)
            test_point = np.array([x, y, z])

            # Find minimum distance to any seed point
            min_dist = float("inf")
            for _, _, px, py, pz in points:
                seed_point = np.array([px, py, pz])
                dist = self.spherical_distance(test_point, seed_point)
                min_dist = min(min_dist, dist)

            # We want to maximize the minimum distance (negative for minimization)
            return -min_dist

        # Try multiple random starting points to find global maximum
        best_result = None
        best_distance = 0

        for _ in range(20):  # Multiple random starts
            # Random starting point
            theta_start = np.random.uniform(0, 2 * np.pi)
            phi_start = np.random.uniform(0, np.pi)

            try:
                result = minimize(
                    objective,
                    [theta_start, phi_start],
                    method="L-BFGS-B",
                    bounds=[(0, 2 * np.pi), (0, np.pi)],
                )

                if result.success and -result.fun > best_distance:
                    best_distance = -result.fun
                    best_result = result
            except:
                continue

        if best_result is not None:
            theta, phi = best_result.x
            x = np.sin(phi) * np.cos(theta)
            y = np.sin(phi) * np.sin(theta)
            z = np.cos(phi)
            center = np.array([x, y, z])

            # Convert radius from radians to more intuitive units
            # Keep in radians for sphere calculations, but also provide degrees
            radius_radians = best_distance
            radius_degrees = np.degrees(radius_radians)

            return center, radius_radians, radius_degrees

        return None, 0, 0

    def update_points_listbox(self):
        """Update the points listbox display"""
        self.points_listbox.delete(0, tk.END)
        for i, (lat, lon, x, y, z) in enumerate(self.seed_points):
            self.points_listbox.insert(tk.END, f"{i+1}: ({lat:.2f}°, {lon:.2f}°)")

    def display_points_only(self):
        """Display only the seed points on the sphere without Voronoi diagram"""
        if not self.seed_points:
            self.draw_empty_sphere()
            return

        self.ax.clear()

        # Draw sphere wireframe
        u = np.linspace(0, 2 * np.pi, 30)
        v = np.linspace(0, np.pi, 30)
        x_sphere = np.outer(np.cos(u), np.sin(v))
        y_sphere = np.outer(np.sin(u), np.sin(v))
        z_sphere = np.outer(np.ones(np.size(u)), np.cos(v))
        self.ax.plot_wireframe(
            x_sphere, y_sphere, z_sphere, alpha=0.1, color="lightgray"
        )

        # Plot seed points
        points = np.array([[x, y, z] for _, _, x, y, z in self.seed_points])
        self.ax.scatter(
            points[:, 0],
            points[:, 1],
            points[:, 2],
            c="red",
            s=120,
            alpha=0.9,
            label="Seed Points",
            edgecolors="black",
        )

        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.set_zlabel("Z")
        self.ax.set_title(f"Seed Points on Sphere ({len(self.seed_points)} points)")
        self.ax.legend()

        # Set equal aspect ratio
        self.ax.set_box_aspect([1, 1, 1])

        self.canvas.draw()

    def draw_circle_on_sphere(self, center, radius, ax):
        """Draw a circle on the sphere surface"""
        try:
            # Normalize the center point
            center = center / np.linalg.norm(center)

            # Create a circle in the plane perpendicular to the center vector
            # First, find two orthogonal vectors in the plane
            if abs(center[2]) < 0.9:
                v1 = np.cross(center, [0, 0, 1])
            else:
                v1 = np.cross(center, [1, 0, 0])
            v1 = v1 / np.linalg.norm(v1)
            v2 = np.cross(center, v1)
            v2 = v2 / np.linalg.norm(v2)

            # Generate points on the circle
            theta = np.linspace(0, 2 * np.pi, 50)
            circle_points = []

            for t in theta:
                # Point on circle in the tangent plane
                circle_point = center + radius * (np.cos(t) * v1 + np.sin(t) * v2)
                # Project back to sphere surface
                circle_point = circle_point / np.linalg.norm(circle_point)
                circle_points.append(circle_point)

            circle_points = np.array(circle_points)

            # Plot the circle
            ax.plot(
                circle_points[:, 0],
                circle_points[:, 1],
                circle_points[:, 2],
                "g-",
                linewidth=3,
                label="Largest Empty Circle",
            )

            # Plot the center
            ax.scatter(
                [center[0]],
                [center[1]],
                [center[2]],
                c="green",
                s=100,
                marker="*",
                label="Circle Center",
            )

            return True
        except Exception as e:
            print(f"Error drawing circle: {e}")
            return False

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
        self.ax.plot_wireframe(
            x_sphere, y_sphere, z_sphere, alpha=0.3, color="lightgray"
        )

        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.set_zlabel("Z")
        self.ax.set_title("Spherical Surface")

        # Set equal aspect ratio
        self.ax.set_box_aspect([1, 1, 1])

        self.canvas.draw()

    def plot_voronoi(self):
        """Plot the Voronoi diagram on the sphere"""
        if len(self.seed_points) < 4:
            messagebox.showerror(
                "Error", "Need at least 4 seed points to generate Voronoi diagram"
            )
            return

        try:
            # Extract Cartesian coordinates
            points = np.array([[x, y, z] for _, _, x, y, z in self.seed_points])

            print(f"Plotting Voronoi with {len(points)} points")  # Debug

            # Generate spherical Voronoi diagram
            sv = SphericalVoronoi(points, radius=1, center=np.array([0, 0, 0]))
            sv.sort_vertices_of_regions()

            print("Voronoi diagram generated successfully")  # Debug

            # Find largest empty circle
            print("Finding largest empty circle...")  # Debug
            empty_center, radius_rad, radius_deg = self.find_largest_empty_circle(
                self.seed_points
            )
            self.largest_empty_circle = empty_center
            self.empty_circle_radius = radius_rad

            print(
                f"Empty circle found: center={empty_center is not None}, radius={radius_rad}"
            )  # Debug

            # Update the label
            if empty_center is not None and radius_rad > 0:
                self.empty_circle_label.config(
                    text=f"Largest Empty Circle Radius: {radius_deg:.2f}° ({radius_rad:.4f} rad)"
                )
            else:
                self.empty_circle_label.config(
                    text="Largest Empty Circle Radius: Calculation failed"
                )

            self.ax.clear()

            # Draw sphere wireframe
            u = np.linspace(0, 2 * np.pi, 30)
            v = np.linspace(0, np.pi, 30)
            x_sphere = np.outer(np.cos(u), np.sin(v))
            y_sphere = np.outer(np.sin(u), np.sin(v))
            z_sphere = np.outer(np.ones(np.size(u)), np.cos(v))
            self.ax.plot_wireframe(
                x_sphere, y_sphere, z_sphere, alpha=0.1, color="lightgray"
            )

            # Plot seed points
            self.ax.scatter(
                points[:, 0],
                points[:, 1],
                points[:, 2],
                c="red",
                s=120,
                alpha=0.9,
                label="Seed Points",
                edgecolors="black",
            )

            # Plot Voronoi vertices
            if hasattr(sv, "vertices") and len(sv.vertices) > 0:
                print(f"Plotting {len(sv.vertices)} Voronoi vertices")  # Debug
                self.ax.scatter(
                    sv.vertices[:, 0],
                    sv.vertices[:, 1],
                    sv.vertices[:, 2],
                    c="blue",
                    s=60,
                    alpha=0.7,
                    label="Voronoi Vertices",
                )

            # Plot Voronoi edges with better error handling
            if hasattr(sv, "regions") and len(sv.regions) > 0:
                print(f"Plotting {len(sv.regions)} Voronoi regions")  # Debug
                colors = plt.cm.Set3(np.linspace(0, 1, len(sv.regions)))

                edge_count = 0
                for region, color in zip(sv.regions, colors):
                    if len(region) > 2:  # Need at least 3 vertices for a valid region
                        try:
                            # Get vertices of this region
                            region_vertices = sv.vertices[region]

                            # Plot the region boundary
                            for i in range(len(region_vertices)):
                                start = region_vertices[i]
                                end = region_vertices[(i + 1) % len(region_vertices)]
                                self.ax.plot(
                                    [start[0], end[0]],
                                    [start[1], end[1]],
                                    [start[2], end[2]],
                                    color=color,
                                    linewidth=2,
                                    alpha=0.8,
                                )
                                edge_count += 1
                        except (IndexError, ValueError) as e:
                            print(f"Skipping invalid region: {e}")
                            continue
                print(f"Drew {edge_count} Voronoi edges")  # Debug

            # Draw largest empty circle
            if empty_center is not None and radius_rad > 0:
                print("Drawing largest empty circle...")  # Debug
                success = self.draw_circle_on_sphere(empty_center, radius_rad, self.ax)
                if success:
                    print("Empty circle drawn successfully")
                else:
                    print("Failed to draw empty circle")

            self.ax.set_xlabel("X")
            self.ax.set_ylabel("Y")
            self.ax.set_zlabel("Z")
            self.ax.set_title(
                f"Spherical Voronoi Diagram ({len(self.seed_points)} points)"
            )
            self.ax.legend()

            # Set equal aspect ratio
            self.ax.set_box_aspect([1, 1, 1])

            self.canvas.draw()
            print("Plot completed successfully")  # Debug

        except Exception as e:
            error_msg = f"Failed to generate Voronoi diagram: {str(e)}"
            messagebox.showerror("Error", error_msg)
            print(f"Debug: Full error: {e}")  # For debugging
            import traceback

            traceback.print_exc()

    def clear_plot(self):
        """Clear the plot and reset"""
        self.seed_points = []
        self.largest_empty_circle = None
        self.empty_circle_radius = 0.0
        self.empty_circle_label.config(
            text="Largest Empty Circle Radius: Not calculated"
        )
        self.update_points_listbox()
        self.draw_empty_sphere()


def main():
    root = tk.Tk()
    app = SphericalVoronoiApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
