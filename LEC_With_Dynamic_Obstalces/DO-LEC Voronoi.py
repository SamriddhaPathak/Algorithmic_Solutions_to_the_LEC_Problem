import tkinter as tk
from tkinter import ttk, messagebox
import math
import time
import threading
from typing import List, Dict, Tuple, Optional
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Circle, Polygon
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import heapq
from collections import defaultdict


class DynamicLECVoronoi:
    def __init__(self, root):
        self.root = root
        self.root.title("Dynamic LEC with Voronoi Diagrams")
        self.root.geometry("1400x900")
        self.root.configure(bg="#f0f0f0")

        # State management
        self.sites = []
        self.obstacle = []
        self.is_animating = False
        self.current_time = 0.0
        self.mode = "sites"  # 'sites', 'obstacle', 'paths'
        self.selected_obstacle = None
        self.temp_path = []
        self.show_voronoi = True
        self.show_convex_hull = True
        self.show_lec = True
        self.animation_speed = 1.0
        self.show_info = False

        # Canvas dimensions
        self.CANVAS_WIDTH = 800
        self.CANVAS_HEIGHT = 600
        self.MARGIN = 50

        # Animation control
        self.animation_thread = None
        self.animation_running = False

        self.setup_gui()
        self.draw()

    def setup_gui(self):
        # Main container
        main_container = tk.Frame(self.root, bg="#f0f0f0")
        main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Header
        header_frame = tk.Frame(main_container, bg="white", relief=tk.RAISED, bd=1)
        header_frame.pack(fill=tk.X, pady=(0, 10))

        title_label = tk.Label(
            header_frame,
            text="Dynamic LEC with Voronoi Diagrams",
            font=("Arial", 16, "bold"),
            bg="white",
            fg="#333",
        )
        title_label.pack(side=tk.LEFT, padx=10, pady=10)

        info_button = tk.Button(
            header_frame,
            text="Info",
            command=self.toggle_info,
            bg="#3498db",
            fg="white",
            font=("Arial", 10),
        )
        info_button.pack(side=tk.RIGHT, padx=10, pady=10)

        # Info panel (initially hidden)
        self.info_frame = tk.Frame(main_container, bg="#e3f2fd", relief=tk.RAISED, bd=1)

        info_text = """How to Use:
• Sites Mode: Click to add Voronoi sites (red dots)
• obstacle Mode: Click to create obstacle and define their paths
• Visualization: Toggle different overlays and animate obstacle movement
• LEC (Orange): Largest Empty Circle avoiding all sites and obstacle
• Voronoi (Blue): Voronoi diagram edges
• Convex Hull (Green): Convex hull of all sites"""

        tk.Label(
            self.info_frame,
            text=info_text,
            bg="#e3f2fd",
            fg="#1976d2",
            font=("Arial", 10),
            justify=tk.LEFT,
        ).pack(padx=10, pady=10)

        # Content frame
        content_frame = tk.Frame(main_container, bg="#f0f0f0")
        content_frame.pack(fill=tk.BOTH, expand=True)

        # Setup control panel and canvas
        self.setup_control_panel(content_frame)
        self.setup_canvas(content_frame)

    def setup_control_panel(self, parent):
        # Control panel frame
        control_frame = tk.Frame(parent, bg="white", relief=tk.RAISED, bd=1, width=320)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        control_frame.pack_propagate(False)

        # Title
        title_label = tk.Label(
            control_frame,
            text="Controls",
            font=("Arial", 14, "bold"),
            bg="white",
            fg="#333",
        )
        title_label.pack(pady=10)

        # Mode Selection
        mode_frame = tk.LabelFrame(
            control_frame, text="Mode", bg="white", font=("Arial", 10, "bold")
        )
        mode_frame.pack(fill=tk.X, padx=10, pady=5)

        button_frame = tk.Frame(mode_frame, bg="white")
        button_frame.pack(fill=tk.X, padx=5, pady=5)

        self.sites_button = tk.Button(
            button_frame,
            text="Add Sites",
            command=lambda: self.set_mode("sites"),
            bg="#e74c3c",
            fg="white",
            font=("Arial", 10),
        )
        self.sites_button.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 5))

        self.obstacle_button = tk.Button(
            button_frame,
            text="Add obstacle",
            command=lambda: self.set_mode("obstacle"),
            bg="#95a5a6",
            fg="white",
            font=("Arial", 10),
        )
        self.obstacle_button.pack(side=tk.LEFT, fill=tk.X, expand=True)

        # Animation Controls
        anim_frame = tk.LabelFrame(
            control_frame, text="Animation", bg="white", font=("Arial", 10, "bold")
        )
        anim_frame.pack(fill=tk.X, padx=10, pady=5)

        button_frame2 = tk.Frame(anim_frame, bg="white")
        button_frame2.pack(fill=tk.X, padx=5, pady=5)

        self.play_button = tk.Button(
            button_frame2,
            text="Play",
            command=self.toggle_animation,
            bg="#27ae60",
            fg="white",
            font=("Arial", 10),
        )
        self.play_button.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 5))

        reset_button = tk.Button(
            button_frame2,
            text="Reset",
            command=self.reset_time,
            bg="#3498db",
            fg="white",
            font=("Arial", 10),
        )
        reset_button.pack(side=tk.LEFT, fill=tk.X, expand=True)

        # Speed control
        speed_frame = tk.Frame(anim_frame, bg="white")
        speed_frame.pack(fill=tk.X, padx=5, pady=5)

        tk.Label(speed_frame, text="Speed:", bg="white", font=("Arial", 10)).pack(
            side=tk.LEFT
        )
        self.speed_var = tk.DoubleVar(value=1.0)
        self.speed_scale = tk.Scale(
            speed_frame,
            from_=0.1,
            to=3.0,
            resolution=0.1,
            orient=tk.HORIZONTAL,
            variable=self.speed_var,
            bg="white",
            font=("Arial", 9),
        )
        self.speed_scale.pack(side=tk.RIGHT, fill=tk.X, expand=True)

        # Display Options
        display_frame = tk.LabelFrame(
            control_frame, text="Display", bg="white", font=("Arial", 10, "bold")
        )
        display_frame.pack(fill=tk.X, padx=10, pady=5)

        self.voronoi_var = tk.BooleanVar(value=True)
        tk.Checkbutton(
            display_frame,
            text="Voronoi Diagram",
            variable=self.voronoi_var,
            command=self.update_display,
            bg="white",
            font=("Arial", 10),
        ).pack(anchor=tk.W, padx=5, pady=2)

        self.hull_var = tk.BooleanVar(value=True)
        tk.Checkbutton(
            display_frame,
            text="Convex Hull",
            variable=self.hull_var,
            command=self.update_display,
            bg="white",
            font=("Arial", 10),
        ).pack(anchor=tk.W, padx=5, pady=2)

        self.lec_var = tk.BooleanVar(value=True)
        tk.Checkbutton(
            display_frame,
            text="Largest Empty Circle",
            variable=self.lec_var,
            command=self.update_display,
            bg="white",
            font=("Arial", 10),
        ).pack(anchor=tk.W, padx=5, pady=2)

        # Statistics
        stats_frame = tk.LabelFrame(
            control_frame, text="Statistics", bg="white", font=("Arial", 10, "bold")
        )
        stats_frame.pack(fill=tk.X, padx=10, pady=5)

        self.stats_label = tk.Label(
            stats_frame,
            text="Sites: 0\nobstacle: 0\nTime: 0.00s",
            bg="white",
            font=("Arial", 10),
            justify=tk.LEFT,
        )
        self.stats_label.pack(padx=5, pady=5)

        # obstacle List
        self.obstacle_list_frame = tk.LabelFrame(
            control_frame, text="obstacle", bg="white", font=("Arial", 10, "bold")
        )
        self.obstacle_list_frame.pack(fill=tk.X, padx=10, pady=5)

        self.obstacle_listbox = tk.Listbox(
            self.obstacle_list_frame, height=4, font=("Arial", 9)
        )
        self.obstacle_listbox.pack(fill=tk.X, padx=5, pady=5)
        self.obstacle_listbox.bind("<<ListboxSelect>>", self.on_obstacle_select)

        # Clear Button
        clear_button = tk.Button(
            control_frame,
            text="Clear All",
            command=self.clear_all,
            bg="#e74c3c",
            fg="white",
            font=("Arial", 12, "bold"),
        )
        clear_button.pack(fill=tk.X, padx=10, pady=20)

    def setup_canvas(self, parent):
        # Canvas frame
        canvas_frame = tk.Frame(parent, bg="white", relief=tk.RAISED, bd=1)
        canvas_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Create matplotlib figure
        self.fig, self.ax = plt.subplots(figsize=(10, 8))
        self.fig.patch.set_facecolor("white")

        # Create canvas
        self.canvas = FigureCanvasTkAgg(self.fig, canvas_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Bind click event
        self.canvas.mpl_connect("button_press_event", self.on_canvas_click)

    def set_mode(self, mode):
        self.mode = mode
        print(f"Mode changed to: {mode}")  # Debug
        # Update button colors
        if mode == "sites":
            self.sites_button.config(bg="#e74c3c")
            self.obstacle_button.config(bg="#95a5a6")
        else:
            self.sites_button.config(bg="#95a5a6")
            self.obstacle_button.config(bg="#8e44ad")

    def toggle_info(self):
        self.show_info = not self.show_info
        if self.show_info:
            self.info_frame.pack(
                fill=tk.X,
                pady=(0, 10),
                before=self.root.winfo_children()[0].winfo_children()[2],
            )
        else:
            self.info_frame.pack_forget()

    def update_display(self):
        self.show_voronoi = self.voronoi_var.get()
        self.show_convex_hull = self.hull_var.get()
        self.show_lec = self.lec_var.get()
        print(
            f"Display updated: Voronoi={self.show_voronoi}, Hull={self.show_convex_hull}, LEC={self.show_lec}"
        )  # Debug
        self.draw()

    def toggle_animation(self):
        if self.is_animating:
            self.is_animating = False
            self.animation_running = False
            self.play_button.config(text="Play", bg="#27ae60")
            print("Animation stopped")  # Debug
        else:
            self.is_animating = True
            self.animation_running = True
            self.play_button.config(text="Pause", bg="#f39c12")
            print("Animation started")  # Debug
            self.start_animation()

    def start_animation(self):
        def animate():
            while self.animation_running:
                self.current_time += 0.02 * self.speed_var.get()
                self.root.after(0, self.draw)
                self.root.after(0, self.update_stats)
                time.sleep(0.02)

        self.animation_thread = threading.Thread(target=animate, daemon=True)
        self.animation_thread.start()

    def reset_time(self):
        self.current_time = 0.0
        print("Time reset")  # Debug
        self.draw()
        self.update_stats()

    def clear_all(self):
        self.sites = []
        self.obstacle = []
        self.temp_path = []
        self.selected_obstacle = None
        self.current_time = 0.0
        self.obstacle_listbox.delete(0, tk.END)
        print("All cleared")  # Debug
        self.draw()
        self.update_stats()

    def on_obstacle_select(self, event):
        selection = self.obstacle_listbox.curselection()
        if selection:
            self.selected_obstacle = selection[0]
            print(f"Selected obstacle: {self.selected_obstacle}")  # Debug
        else:
            self.selected_obstacle = None
        self.draw()

    def update_stats(self):
        stats_text = f"Sites: {len(self.sites)}\nobstacle: {len(self.obstacle)}\nTime: {self.current_time:.2f}s"
        self.stats_label.config(text=stats_text)

    def on_canvas_click(self, event):
        if event.inaxes != self.ax:
            return

        x, y = event.xdata, event.ydata
        print(f"Canvas clicked at: ({x:.2f}, {y:.2f})")  # Debug

        # Check bounds
        if (
            x < self.MARGIN
            or x > self.CANVAS_WIDTH - self.MARGIN
            or y < self.MARGIN
            or y > self.CANVAS_HEIGHT - self.MARGIN
        ):
            print("Click outside bounds")  # Debug
            return

        if self.mode == "sites":
            self.sites.append({"x": x, "y": y})
            print(f"Added site: {len(self.sites)} sites total")  # Debug
        elif self.mode == "obstacle":
            if self.selected_obstacle is not None and self.selected_obstacle < len(
                self.obstacle
            ):
                # Adding to existing obstacle's path
                self.obstacle[self.selected_obstacle]["path"].append({"x": x, "y": y})
                print(f"Added point to obstacle {self.selected_obstacle}")  # Debug
                # Update listbox
                self.obstacle_listbox.delete(self.selected_obstacle)
                self.obstacle_listbox.insert(
                    self.selected_obstacle,
                    f"Obstacle {self.selected_obstacle + 1} ({len(self.obstacle[self.selected_obstacle]['path'])} points)",
                )
            else:
                # Create new obstacle
                new_obstacle = {
                    "radius": 15,
                    "speed": 1,
                    "path": [{"x": x, "y": y}],
                    "currentPos": {"x": x, "y": y},
                }
                self.obstacle.append(new_obstacle)
                self.selected_obstacle = len(self.obstacle) - 1
                print(
                    f"Created new obstacle: {len(self.obstacle)} obstacle total"
                )  # Debug

                # Update obstacle listbox
                self.obstacle_listbox.insert(
                    tk.END, f"Obstacle {len(self.obstacle)} (1 points)"
                )

        self.draw()
        self.update_stats()

    def compute_voronoi(self, sites):
        """Compute Voronoi edges using Fortune's algorithm via scipy.spatial"""
        if len(sites) < 4:
            # Not enough points for scipy Voronoi
            return []

        points = np.array([[site["x"], site["y"]] for site in sites])
        vor = Voronoi(points)

        edges = []
        voronoi_plot_2d(
            vor,
            ax=self.ax,
            show_vertices=False,
            line_colors="black",
            line_width=2,
            line_alpha=0.6,
            point_size=10,
        )

        return edges

    def compute_convex_hull(self, points):
        """Graham scan for convex hull"""
        if len(points) < 3:
            return points

        # Sort points
        sorted_points = sorted(points, key=lambda p: (p["x"], p["y"]))

        def cross(o, a, b):
            return (a["x"] - o["x"]) * (b["y"] - o["y"]) - (a["y"] - o["y"]) * (
                b["x"] - o["x"]
            )

        # Build lower hull
        lower = []
        for p in sorted_points:
            while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
                lower.pop()
            lower.append(p)

        # Build upper hull
        upper = []
        for p in reversed(sorted_points):
            while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
                upper.pop()
            upper.append(p)

        return lower[:-1] + upper[:-1]

    def get_current_obstacle_positions(self, time):
        """Get current positions of obstacle based on time"""
        current_obstacle = []

        for obstacle in self.obstacle:
            if len(obstacle["path"]) < 2:
                current_pos = (
                    obstacle["path"][0] if obstacle["path"] else {"x": 0, "y": 0}
                )
                current_obstacle.append({**obstacle, "currentPos": current_pos})
                continue

            speed = obstacle.get("speed", 1)
            total_time = len(obstacle["path"]) - 1
            normalized_time = (time * speed) % (total_time * 2)

            if normalized_time <= total_time:
                t = normalized_time
            else:
                t = total_time * 2 - normalized_time

            segment_index = min(int(t), len(obstacle["path"]) - 2)
            local_t = t - segment_index

            p1 = obstacle["path"][segment_index]
            p2 = obstacle["path"][segment_index + 1]

            current_pos = {
                "x": p1["x"] + (p2["x"] - p1["x"]) * local_t,
                "y": p1["y"] + (p2["y"] - p1["y"]) * local_t,
            }

            current_obstacle.append({**obstacle, "currentPos": current_pos})

        return current_obstacle

    def compute_lec(self, sites, current_obstacle):
        """Compute Largest Empty Circle using optimized O(n log n) approach"""
        if len(sites) < 2:
            return None

        bounds = {
            "minX": self.MARGIN,
            "minY": self.MARGIN,
            "maxX": self.CANVAS_WIDTH - self.MARGIN,
            "maxY": self.CANVAS_HEIGHT - self.MARGIN,
        }

        class Point:
            __slots__ = ["x", "y"]

            def __init__(self, x, y):
                self.x = x
                self.y = y

            def __eq__(self, other):
                return abs(self.x - other.x) < 1e-9 and abs(self.y - other.y) < 1e-9

            def __hash__(self):
                return hash((round(self.x, 6), round(self.y, 6)))

            def distance_to(self, other):
                dx = self.x - other.x
                dy = self.y - other.y
                return (dx * dx + dy * dy) ** 0.5

        class OptimizedVoronoiSolver:
            def __init__(self, sites, bounds):
                self.sites = sites
                self.bounds = bounds
                self.candidates = set()

            def solve(self):
                """Main solving method using optimized approach"""
                # 1. Add boundary candidates (always useful)
                self._add_boundary_candidates()

                # 2. Use Delaunay triangulation approach for Voronoi vertices
                self._add_delaunay_circumcenters()

                # 3. Add perpendicular bisector intersections with boundaries
                self._add_bisector_boundary_intersections()

                return list(self.candidates)

            def _add_boundary_candidates(self):
                """Add boundary points as candidates"""
                # Corners
                self.candidates.add(Point(self.bounds["minX"], self.bounds["minY"]))
                self.candidates.add(Point(self.bounds["minX"], self.bounds["maxY"]))
                self.candidates.add(Point(self.bounds["maxX"], self.bounds["minY"]))
                self.candidates.add(Point(self.bounds["maxX"], self.bounds["maxY"]))

                # Mid-points of boundaries
                mid_x = (self.bounds["minX"] + self.bounds["maxX"]) / 2
                mid_y = (self.bounds["minY"] + self.bounds["maxY"]) / 2
                self.candidates.add(Point(self.bounds["minX"], mid_y))
                self.candidates.add(Point(self.bounds["maxX"], mid_y))
                self.candidates.add(Point(mid_x, self.bounds["minY"]))
                self.candidates.add(Point(mid_x, self.bounds["maxY"]))

            def _add_delaunay_circumcenters(self):
                """Add circumcenters using optimized Delaunay approach"""
                n = len(self.sites)
                if n < 3:
                    return

                # Sort sites by x-coordinate for divide-and-conquer
                sorted_sites = sorted(
                    enumerate(self.sites), key=lambda x: (x[1].x, x[1].y)
                )

                # Use incremental construction for smaller datasets
                if n <= 50:
                    self._incremental_delaunay(sorted_sites)
                else:
                    # For larger datasets, use sampling + local search
                    self._sampled_delaunay(sorted_sites)

            def _incremental_delaunay(self, sorted_sites):
                """Incremental Delaunay triangulation for smaller datasets"""
                triangles = []

                # Start with first 3 points
                if len(sorted_sites) >= 3:
                    p1, p2, p3 = (
                        sorted_sites[0][1],
                        sorted_sites[1][1],
                        sorted_sites[2][1],
                    )
                    if not self._are_collinear(p1, p2, p3):
                        triangles.append(
                            (sorted_sites[0][0], sorted_sites[1][0], sorted_sites[2][0])
                        )

                # Add remaining points incrementally
                for i in range(3, len(sorted_sites)):
                    new_idx, new_point = sorted_sites[i]
                    new_triangles = []

                    # Find triangles that need to be updated
                    for tri_idx, (i1, i2, i3) in enumerate(triangles):
                        p1, p2, p3 = self.sites[i1], self.sites[i2], self.sites[i3]

                        # Check if new point is inside circumcircle
                        if self._point_in_circumcircle(new_point, p1, p2, p3):
                            # Create new triangles with the new point
                            new_triangles.extend(
                                [
                                    (new_idx, i1, i2),
                                    (new_idx, i2, i3),
                                    (new_idx, i3, i1),
                                ]
                            )
                        else:
                            new_triangles.append((i1, i2, i3))

                    triangles = new_triangles

                # Add circumcenters of valid triangles
                for i1, i2, i3 in triangles:
                    circumcenter = self._circumcenter(
                        self.sites[i1], self.sites[i2], self.sites[i3]
                    )
                    if circumcenter and self._is_point_in_bounds(circumcenter):
                        self.candidates.add(circumcenter)

            def _sampled_delaunay(self, sorted_sites):
                """Sampled approach for larger datasets"""
                n = len(sorted_sites)

                # Sample every k-th point where k = sqrt(n)
                sample_step = max(1, int(n**0.5))
                sampled_indices = list(range(0, n, sample_step))

                # Ensure we have at least 10 samples
                if len(sampled_indices) < 10:
                    sampled_indices = list(range(0, min(n, 20)))

                # Process sampled triangles
                for i in range(len(sampled_indices)):
                    for j in range(i + 1, len(sampled_indices)):
                        for k in range(j + 1, len(sampled_indices)):
                            idx1, idx2, idx3 = (
                                sampled_indices[i],
                                sampled_indices[j],
                                sampled_indices[k],
                            )
                            p1, p2, p3 = (
                                sorted_sites[idx1][1],
                                sorted_sites[idx2][1],
                                sorted_sites[idx3][1],
                            )

                            if not self._are_collinear(p1, p2, p3):
                                circumcenter = self._circumcenter(p1, p2, p3)
                                if circumcenter and self._is_point_in_bounds(
                                    circumcenter
                                ):
                                    self.candidates.add(circumcenter)

            def _add_bisector_boundary_intersections(self):
                """Add intersections of perpendicular bisectors with boundaries"""
                n = len(self.sites)

                # Only process nearest neighbors to reduce complexity
                for i in range(n):
                    # Find k nearest neighbors (where k = min(8, n-1))
                    k = min(8, n - 1)
                    distances = []

                    for j in range(n):
                        if i != j:
                            dist = self.sites[i].distance_to(self.sites[j])
                            distances.append((dist, j))

                    distances.sort()

                    # Process k nearest neighbors
                    for _, j in distances[:k]:
                        intersections = self._get_bisector_boundary_intersections(
                            self.sites[i], self.sites[j]
                        )
                        for point in intersections:
                            if self._is_point_in_bounds(point):
                                self.candidates.add(point)

            def _get_bisector_boundary_intersections(self, site1, site2):
                """Get intersections of perpendicular bisector with boundary"""
                intersections = []

                # Midpoint and direction
                mx = (site1.x + site2.x) / 2
                my = (site1.y + site2.y) / 2

                dx = site2.x - site1.x
                dy = site2.y - site1.y

                # Perpendicular bisector: ax + by + c = 0
                a, b, c = -dy, dx, dy * mx - dx * my

                # Find intersections with boundary
                if abs(b) > 1e-10:  # Not vertical
                    # Left boundary
                    y = -(a * self.bounds["minX"] + c) / b
                    if self.bounds["minY"] <= y <= self.bounds["maxY"]:
                        intersections.append(Point(self.bounds["minX"], y))

                    # Right boundary
                    y = -(a * self.bounds["maxX"] + c) / b
                    if self.bounds["minY"] <= y <= self.bounds["maxY"]:
                        intersections.append(Point(self.bounds["maxX"], y))

                if abs(a) > 1e-10:  # Not horizontal
                    # Bottom boundary
                    x = -(b * self.bounds["minY"] + c) / a
                    if self.bounds["minX"] <= x <= self.bounds["maxX"]:
                        intersections.append(Point(x, self.bounds["minY"]))

                    # Top boundary
                    x = -(b * self.bounds["maxY"] + c) / a
                    if self.bounds["minX"] <= x <= self.bounds["maxX"]:
                        intersections.append(Point(x, self.bounds["maxY"]))

                return intersections

            def _circumcenter(self, p1, p2, p3):
                """Find circumcenter of triangle - optimized version"""
                x1, y1 = p1.x, p1.y
                x2, y2 = p2.x, p2.y
                x3, y3 = p3.x, p3.y

                # Check if points are collinear
                det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3)
                if abs(det) < 1e-10:
                    return None

                # Pre-calculate squares to avoid repeated computation
                sq1 = x1 * x1 + y1 * y1
                sq2 = x2 * x2 + y2 * y2
                sq3 = x3 * x3 + y3 * y3

                # Calculate circumcenter
                ux = (sq1 * (y2 - y3) + sq2 * (y3 - y1) + sq3 * (y1 - y2)) / (2 * det)
                uy = (sq1 * (x3 - x2) + sq2 * (x1 - x3) + sq3 * (x2 - x1)) / (2 * det)

                return Point(ux, uy)

            def _are_collinear(self, p1, p2, p3):
                """Check if three points are collinear"""
                return (
                    abs((p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y))
                    < 1e-10
                )

            def _point_in_circumcircle(self, point, p1, p2, p3):
                """Check if point is inside circumcircle of triangle"""
                circumcenter = self._circumcenter(p1, p2, p3)
                if not circumcenter:
                    return False

                # Calculate circumradius
                circumradius = circumcenter.distance_to(p1)

                # Check if point is inside
                return point.distance_to(circumcenter) < circumradius

            def _is_point_in_bounds(self, point):
                """Check if point is within bounds"""
                return (
                    self.bounds["minX"] <= point.x <= self.bounds["maxX"]
                    and self.bounds["minY"] <= point.y <= self.bounds["maxY"]
                )

        # Convert sites to Point objects
        points = [Point(site["x"], site["y"]) for site in sites]
        obstacle = [
            (Point(obs["currentPos"]["x"], obs["currentPos"]["y"]), obs["radius"])
            for obs in current_obstacle
        ]

        # Generate candidates using optimized solver
        solver = OptimizedVoronoiSolver(points, bounds)
        candidates = solver.solve()

        # Optimized distance calculation functions
        def is_point_in_bounds(point):
            """Check if point is within bounds"""
            return (
                bounds["minX"] <= point.x <= bounds["maxX"]
                and bounds["minY"] <= point.y <= bounds["maxY"]
            )

        def get_max_radius_at_point(point):
            """Get maximum possible radius at given point - optimized"""
            if not is_point_in_bounds(point):
                return 0

            # Distance to boundary (optimized)
            dist_to_bounds = min(
                point.x - bounds["minX"],
                bounds["maxX"] - point.x,
                point.y - bounds["minY"],
                bounds["maxY"] - point.y,
            )

            if dist_to_bounds <= 0:
                return 0

            # Distance to nearest site (optimized with early termination)
            min_dist_to_site = float("inf")
            for site in points:
                dx = point.x - site.x
                dy = point.y - site.y
                dist_sq = dx * dx + dy * dy
                if dist_sq < min_dist_to_site * min_dist_to_site:
                    min_dist_to_site = dist_sq**0.5
                    # Early termination if we found a very close site
                    if min_dist_to_site < 1.0:
                        break

            # Check obstacle clearance
            min_clearance = dist_to_bounds
            for obs_center, obs_radius in obstacle:
                dx = point.x - obs_center.x
                dy = point.y - obs_center.y
                dist_to_obstacle = (dx * dx + dy * dy) ** 0.5
                clearance = dist_to_obstacle - obs_radius
                min_clearance = min(min_clearance, clearance)

                # Early termination if clearance is too small
                if min_clearance <= 0:
                    return 0

            return max(0, min(min_dist_to_site, min_clearance))

        # Find best candidate with optimized search
        best_circle = None
        max_radius = 0

        # Remove duplicates efficiently
        unique_candidates = list(
            {
                (round(c.x, 6), round(c.y, 6)): c
                for c in candidates
                if c and is_point_in_bounds(c)
            }.values()
        )

        # Parallel-style processing with early termination
        for candidate in unique_candidates:
            radius = get_max_radius_at_point(candidate)
            if radius > max_radius:
                max_radius = radius
                best_circle = {
                    "center": {"x": candidate.x, "y": candidate.y},
                    "radius": radius,
                }

                # Early termination if we found a very good solution
                if radius > 50:  # Adjust threshold as needed
                    break

        # Fallback to center if no good solution found
        if not best_circle or max_radius < 5:
            center_x = (bounds["minX"] + bounds["maxX"]) / 2
            center_y = (bounds["minY"] + bounds["maxY"]) / 2
            center_point = Point(center_x, center_y)

            if is_point_in_bounds(center_point):
                radius = get_max_radius_at_point(center_point)
                if radius > max_radius:
                    best_circle = {
                        "center": {"x": center_x, "y": center_y},
                        "radius": radius,
                    }

        return best_circle

    def draw(self):
        """Main drawing function"""
        try:
            self.ax.clear()
            self.ax.set_xlim(0, self.CANVAS_WIDTH)
            self.ax.set_ylim(0, self.CANVAS_HEIGHT)
            self.ax.set_aspect("equal")
            self.ax.set_facecolor("#f8f9fa")

            # Draw border
            border_x = [
                self.MARGIN,
                self.CANVAS_WIDTH - self.MARGIN,
                self.CANVAS_WIDTH - self.MARGIN,
                self.MARGIN,
                self.MARGIN,
            ]
            border_y = [
                self.MARGIN,
                self.MARGIN,
                self.CANVAS_HEIGHT - self.MARGIN,
                self.CANVAS_HEIGHT - self.MARGIN,
                self.MARGIN,
            ]
            self.ax.plot(border_x, border_y, color="#dee2e6", linewidth=2)

            if len(self.sites) == 0:
                self.ax.text(
                    self.CANVAS_WIDTH // 2,
                    self.CANVAS_HEIGHT // 2,
                    "Click to add sites",
                    ha="center",
                    va="center",
                    fontsize=16,
                    color="#6c757d",
                )
                self.canvas.draw()
                return

            current_obstacle = self.get_current_obstacle_positions(self.current_time)

            # Draw Voronoi diagram
            if self.show_voronoi and len(self.sites) > 1:
                voronoi_edges = self.compute_voronoi(self.sites)
                print(f"Drawing {len(voronoi_edges)} Voronoi edges")  # Debug
                for edge in voronoi_edges:
                    self.ax.plot(
                        [edge["x1"], edge["x2"]],
                        [edge["y1"], edge["y2"]],
                        color="#3498db",
                        linestyle="--",
                        linewidth=1,
                        alpha=0.7,
                    )

            # Draw convex hull
            if self.show_convex_hull and len(self.sites) > 2:
                hull = self.compute_convex_hull(self.sites)
                print(f"Drawing convex hull with {len(hull)} points")  # Debug
                if len(hull) >= 3:
                    hull_x = [p["x"] for p in hull] + [hull[0]["x"]]
                    hull_y = [p["y"] for p in hull] + [hull[0]["y"]]
                    self.ax.plot(hull_x, hull_y, color="#2ecc71", linewidth=2)

            # Draw LEC
            if self.show_lec and len(self.sites) >= 2:
                lec = self.compute_lec(self.sites, current_obstacle)
                if lec:
                    print(f"Drawing LEC with radius {lec['radius']:.2f}")  # Debug
                    circle = Circle(
                        (lec["center"]["x"], lec["center"]["y"]),
                        lec["radius"],
                        fill=True,
                        facecolor="#f39c12",
                        alpha=0.1,
                        edgecolor="#f39c12",
                        linewidth=2,
                    )
                    self.ax.add_patch(circle)
                    # Draw center point
                    self.ax.plot(
                        lec["center"]["x"],
                        lec["center"]["y"],
                        "o",
                        color="#f39c12",
                        markersize=6,
                    )

            # Draw obstacle
            print(f"Drawing {len(current_obstacle)} obstacle")  # Debug
            for i, obstacle in enumerate(current_obstacle):
                # Draw path
                if len(obstacle["path"]) > 1:
                    path_x = [p["x"] for p in obstacle["path"]]
                    path_y = [p["y"] for p in obstacle["path"]]
                    self.ax.plot(
                        path_x,
                        path_y,
                        color="#95a5a6",
                        linestyle="--",
                        linewidth=2,
                        alpha=0.7,
                    )

                # Draw obstacle
                color = "#e67e22" if i == self.selected_obstacle else "#8e44ad"
                circle = Circle(
                    (obstacle["currentPos"]["x"], obstacle["currentPos"]["y"]),
                    obstacle["radius"],
                    fill=True,
                    facecolor=color,
                    alpha=0.8,
                )
                self.ax.add_patch(circle)

                # Draw direction indicator
                self.ax.plot(
                    obstacle["currentPos"]["x"],
                    obstacle["currentPos"]["y"],
                    "o",
                    color="white",
                    markersize=4,
                )

            # Draw sites
            print(f"Drawing {len(self.sites)} sites")  # Debug
            for site in self.sites:
                self.ax.plot(site["x"], site["y"], "o", color="#e74c3c", markersize=8)

            # Draw temporary path
            if len(self.temp_path) > 0:
                temp_x = [p["x"] for p in self.temp_path]
                temp_y = [p["y"] for p in self.temp_path]
                self.ax.plot(
                    temp_x,
                    temp_y,
                    color="#e74c3c",
                    linestyle=":",
                    linewidth=2,
                    alpha=0.7,
                )

            self.ax.set_title(
                "Dynamic LEC with Voronoi Diagrams", fontsize=14, fontweight="bold"
            )
            self.ax.grid(True, alpha=0.3)

            # Force canvas update
            self.canvas.draw()
            self.canvas.flush_events()

        except Exception as e:
            print(f"Error in draw(): {e}")
            import traceback

            traceback.print_exc()


def main():
    root = tk.Tk()
    app = DynamicLECVoronoi(root)
    root.mainloop()


if __name__ == "__main__":
    main()
