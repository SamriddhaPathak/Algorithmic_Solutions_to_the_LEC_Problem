import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from scipy.spatial import Voronoi, ConvexHull, distance
from scipy.optimize import minimize
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
import warnings
warnings.filterwarnings('ignore')

class InteractiveVoronoiPlotter:
    """
    Interactive Voronoi Diagram Plotter - Start with empty graph
    
    Controls:
    - Left click: Add site point
    - Right click: Add obstacle
    - Middle click or Ctrl+Left click: Remove nearest point
    - Press 'c': Clear all points
    - Press 'r': Generate random points
    - Press 'v': Toggle Voronoi diagram
    - Press 'h': Toggle convex hull
    - Press 'e': Toggle largest empty circle
    - Press 's': Save plot
    """
    
    def __init__(self, width=10, height=8, figsize=(12, 10)):
        self.width = width
        self.height = height
        self.figsize = figsize
        
        # Data storage - start empty
        self.sites = []
        self.obstacles = []
        self.voronoi = None
        self.convex_hull = None
        self.largest_empty_circle = None
        
        # Visualization settings
        self.show_voronoi = True
        self.show_convex_hull = True
        self.show_largest_circle = True
        self.show_obstacles = True
        self.show_sites = True
        
        # Colors
        self.colors = {
            'sites': '#2E4057',
            'obstacles': '#E74C3C',
            'voronoi': '#3498DB',
            'convex_hull': '#E74C3C',
            'largest_circle': '#27AE60',
            'largest_circle_fill': '#27AE60'
        }
        
        # Initialize plot
        self.fig, self.ax = plt.subplots(figsize=self.figsize)
        self.setup_plot()
        self.setup_event_handlers()
        
        # Show initial empty plot
        self.update_plot()
    
    def setup_plot(self):
        """Setup the matplotlib plot with styling"""
        self.ax.set_xlim(0, self.width)
        self.ax.set_ylim(0, self.height)
        self.ax.set_aspect('equal')
        self.ax.grid(True, alpha=0.3)
        self.ax.set_facecolor('#F8F9FA')
        
        # Styling
        self.ax.set_title('Interactive Voronoi Diagram Plotter\n'
                         'Left Click: Site | Right Click: Obstacle | Press H for Help',
                         fontsize=16, fontweight='bold', pad=20)
        self.ax.set_xlabel('X Coordinate', fontsize=12)
        self.ax.set_ylabel('Y Coordinate', fontsize=12)
        
        plt.tight_layout()
    
    def setup_event_handlers(self):
        """Setup event handlers for interactivity"""
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        
        # Print instructions
        print("\n" + "="*60)
        print("INTERACTIVE VORONOI DIAGRAM PLOTTER")
        print("="*60)
        print("Mouse Controls:")
        print("  • Left Click:           Add site point (blue)")
        print("  • Right Click:          Add obstacle (red X)")
        print("  • Middle Click:         Remove nearest point")
        print("  • Ctrl+Left Click:      Remove nearest point")
        print("\nKeyboard Controls:")
        print("  • 'c' or 'C':          Clear all points")
        print("  • 'r' or 'R':          Generate random points")
        print("  • 'v' or 'V':          Toggle Voronoi diagram")
        print("  • 'h' or 'H':          Toggle convex hull")
        print("  • 'e' or 'E':          Toggle largest empty circle")
        print("  • 's' or 'S':          Save current plot")
        print("  • 'help':              Show this help again")
        print("\nAlgorithms:")
        print("  • Voronoi Diagram:     Fortune's Algorithm O(n log n)")
        print("  • Convex Hull:         Graham Scan O(n log n)")
        print("  • Largest Empty Circle: Optimization approach")
        print("="*60)
        print("Start by clicking to add points!")
        print("="*60 + "\n")
    
    def add_site(self, x, y):
        """Add a site point"""
        if 0 <= x <= self.width and 0 <= y <= self.height:
            self.sites.append([x, y])
            print(f"Added site at ({x:.2f}, {y:.2f})")
            self.compute_all()
            self.update_plot()
    
    def add_obstacle(self, x, y):
        """Add an obstacle point"""
        if 0 <= x <= self.width and 0 <= y <= self.height:
            self.obstacles.append([x, y])
            print(f"Added obstacle at ({x:.2f}, {y:.2f})")
            self.compute_all()
            self.update_plot()
    
    def remove_nearest_point(self, x, y):
        """Remove the nearest point (site or obstacle) to the click location"""
        if not self.sites and not self.obstacles:
            return
            
        all_points = []
        point_types = []
        
        # Collect all points with their types
        for site in self.sites:
            all_points.append(site)
            point_types.append('site')
        for obstacle in self.obstacles:
            all_points.append(obstacle)
            point_types.append('obstacle')
        
        if not all_points:
            return
        
        # Find nearest point
        click_point = np.array([x, y])
        distances = [np.linalg.norm(np.array(point) - click_point) for point in all_points]
        nearest_idx = np.argmin(distances)
        
        # Remove the nearest point
        if point_types[nearest_idx] == 'site':
            site_idx = sum(1 for i in range(nearest_idx) if point_types[i] == 'site')
            removed_point = self.sites.pop(site_idx)
            print(f"Removed site at ({removed_point[0]:.2f}, {removed_point[1]:.2f})")
        else:
            obstacle_idx = sum(1 for i in range(nearest_idx) if point_types[i] == 'obstacle')
            removed_point = self.obstacles.pop(obstacle_idx)
            print(f"Removed obstacle at ({removed_point[0]:.2f}, {removed_point[1]:.2f})")
        
        self.compute_all()
        self.update_plot()
    
    def clear_all(self):
        """Clear all points and obstacles"""
        self.sites = []
        self.obstacles = []
        self.voronoi = None
        self.convex_hull = None
        self.largest_empty_circle = None
        print("Cleared all points")
        self.update_plot()
    
    def generate_random_points(self, n_sites=8, n_obstacles=3):
        """Generate random sites and obstacles"""
        self.clear_all()
        
        # Generate random sites with some spacing to avoid numerical issues
        np.random.seed(None)  # Use random seed each time
        margin = 0.8
        sites = []
        
        # Generate well-spaced sites
        for _ in range(n_sites):
            attempts = 0
            while attempts < 50:  # Prevent infinite loops
                x = np.random.uniform(margin, self.width - margin)
                y = np.random.uniform(margin, self.height - margin)
                new_site = [x, y]
                
                # Check minimum distance from existing sites
                too_close = False
                for existing_site in sites:
                    if np.linalg.norm(np.array(new_site) - np.array(existing_site)) < 0.5:
                        too_close = True
                        break
                
                if not too_close:
                    sites.append(new_site)
                    break
                attempts += 1
            
            if attempts >= 50:  # If we couldn't place it, just add it anyway
                sites.append([np.random.uniform(margin, self.width - margin),
                             np.random.uniform(margin, self.height - margin)])
        
        self.sites = sites
        
        # Generate random obstacles
        obstacles = []
        for _ in range(n_obstacles):
            x = np.random.uniform(margin, self.width - margin)
            y = np.random.uniform(margin, self.height - margin)
            obstacles.append([x, y])
        
        self.obstacles = obstacles
        
        print(f"Generated {len(self.sites)} random sites and {len(self.obstacles)} obstacles")
        self.compute_all()
        self.update_plot()
    
    def compute_voronoi(self):
        """Compute Voronoi diagram using Fortune's Algorithm (scipy implementation)"""
        if len(self.sites) < 2:
            self.voronoi = None
            return
        
        try:
            sites_array = np.array(self.sites)
            # Add small random noise to avoid degenerate cases
            if len(sites_array) > 1:
                # Check for duplicate points
                unique_sites = []
                for site in sites_array:
                    is_duplicate = False
                    for existing in unique_sites:
                        if np.linalg.norm(site - existing) < 1e-10:
                            is_duplicate = True
                            break
                    if not is_duplicate:
                        unique_sites.append(site)
                
                if len(unique_sites) >= 2:
                    sites_array = np.array(unique_sites)
                    self.voronoi = Voronoi(sites_array)
                else:
                    self.voronoi = None
            else:
                self.voronoi = None
        except Exception as e:
            print(f"Error computing Voronoi diagram: {e}")
            self.voronoi = None
    
    def compute_convex_hull(self):
        """Compute convex hull using Graham Scan (scipy implementation)"""
        if len(self.sites) < 3:
            self.convex_hull = None
            return
        
        try:
            sites_array = np.array(self.sites)
            hull = ConvexHull(sites_array)
            self.convex_hull = hull
        except Exception as e:
            print(f"Error computing convex hull: {e}")
            self.convex_hull = None
    
    def compute_largest_empty_circle(self):
        """
        Compute the largest empty circle that doesn't intersect with sites or obstacles
        Uses optimization approach for accuracy
        """
        if len(self.sites) < 2:
            self.largest_empty_circle = None
            return
        
        try:
            # Combine sites and obstacles
            all_points = np.array(self.sites + self.obstacles)
            
            # Define the objective function (negative radius to maximize)
            def objective(center):
                distances = distance.cdist([center], all_points)[0]
                return -np.min(distances)  # Negative for maximization
            
            # Define bounds (within the plot area)
            bounds = [(0.1, self.width-0.1), (0.1, self.height-0.1)]
            
            # Try multiple starting points for global optimization
            best_result = None
            best_radius = 0
            
            # Grid search for initial points
            for start_x in np.linspace(1, self.width-1, 5):
                for start_y in np.linspace(1, self.height-1, 5):
                    try:
                        result = minimize(objective, [start_x, start_y], 
                                        bounds=bounds, method='L-BFGS-B')
                        
                        if result.success:
                            radius = -result.fun
                            if radius > best_radius:
                                best_radius = radius
                                best_result = result
                    except:
                        continue
            
            if best_result and best_radius > 0.01:
                self.largest_empty_circle = {
                    'center': best_result.x,
                    'radius': best_radius
                }
            else:
                self.largest_empty_circle = None
                
        except Exception as e:
            print(f"Error computing largest empty circle: {e}")
            self.largest_empty_circle = None
    
    def compute_all(self):
        """Compute all geometric structures"""
        self.compute_voronoi()
        self.compute_convex_hull()
        self.compute_largest_empty_circle()
    
    def plot_voronoi(self):
        """Plot Voronoi diagram"""
        if not self.show_voronoi or self.voronoi is None:
            return
        
        try:
            # Plot finite Voronoi edges
            for simplex in self.voronoi.ridge_vertices:
                simplex = np.asarray(simplex)  # Ensure it's a numpy array
                if len(simplex) == 2 and -1 not in simplex:  # Finite edge
                    points = self.voronoi.vertices[simplex]
                    # Check if points are within reasonable bounds
                    if np.all(np.abs(points) < 1000):
                        self.ax.plot(points[:, 0], points[:, 1], 
                                   color=self.colors['voronoi'], linewidth=1.5, alpha=0.8)
            
            # Plot infinite edges (clipped to plot bounds) - simplified approach
            center = np.mean(self.voronoi.points, axis=0)
            for pointidx, simplex in zip(self.voronoi.ridge_points, self.voronoi.ridge_vertices):
                simplex = np.asarray(simplex)  # Ensure it's a numpy array
                if -1 in simplex and len(simplex) == 2:  # Infinite edge
                    try:
                        finite_vertices = simplex[simplex >= 0]
                        if len(finite_vertices) > 0:
                            i = finite_vertices[0]  # Finite vertex
                            
                            # Direction vector for the infinite edge
                            t = self.voronoi.points[pointidx[1]] - self.voronoi.points[pointidx[0]]
                            t_norm = np.linalg.norm(t)
                            if t_norm > 0:
                                t = t / t_norm
                                n = np.array([-t[1], t[0]])  # Perpendicular vector
                                
                                # Calculate far point within plot bounds
                                vertex = self.voronoi.vertices[i]
                                if np.all(np.abs(vertex) < 1000):  # Sanity check
                                    # Extend line to plot boundary
                                    extend_distance = max(self.width, self.height) * 2
                                    far_point = vertex + n * extend_distance
                                    
                                    # Clip to plot bounds
                                    if (0 <= vertex[0] <= self.width and 0 <= vertex[1] <= self.height):
                                        self.ax.plot([vertex[0], far_point[0]], 
                                                   [vertex[1], far_point[1]], 
                                                   color=self.colors['voronoi'], 
                                                   linewidth=1.5, alpha=0.6)
                    except (IndexError, ValueError):
                        continue  # Skip problematic edges
                        
        except Exception as e:
            print(f"Warning: Could not plot some Voronoi edges: {e}")
            # Fallback: just plot finite edges
            try:
                for simplex in self.voronoi.ridge_vertices:
                    simplex = np.asarray(simplex)
                    if len(simplex) == 2 and -1 not in simplex:
                        points = self.voronoi.vertices[simplex]
                        if np.all(np.isfinite(points)) and np.all(np.abs(points) < 1000):
                            self.ax.plot(points[:, 0], points[:, 1], 
                                       color=self.colors['voronoi'], linewidth=1.5, alpha=0.8)
            except:
                pass  # If all else fails, skip Voronoi plotting
    
    def plot_convex_hull(self):
        """Plot convex hull"""
        if not self.show_convex_hull or self.convex_hull is None:
            return
        
        # Get hull vertices
        hull_points = np.array(self.sites)[self.convex_hull.vertices]
        
        # Close the hull
        hull_points = np.vstack([hull_points, hull_points[0]])
        
        self.ax.plot(hull_points[:, 0], hull_points[:, 1], 
                    color=self.colors['convex_hull'], linewidth=3, 
                    linestyle='--', alpha=0.9, label='Convex Hull')
    
    def plot_largest_empty_circle(self):
        """Plot largest empty circle"""
        if not self.show_largest_circle or self.largest_empty_circle is None:
            return
        
        center = self.largest_empty_circle['center']
        radius = self.largest_empty_circle['radius']
        
        # Plot circle
        circle = Circle(center, radius, fill=True, 
                       facecolor=self.colors['largest_circle_fill'], alpha=0.2,
                       edgecolor=self.colors['largest_circle'], linewidth=2,
                       linestyle=':', label=f'Largest Empty Circle (r={radius:.2f})')
        self.ax.add_patch(circle)
        
        # Plot center
        self.ax.plot(center[0], center[1], 'o', 
                    color=self.colors['largest_circle'], markersize=8, 
                    markeredgecolor='white', markeredgewidth=2)
    
    def plot_sites(self):
        """Plot site points"""
        if not self.show_sites or not self.sites:
            return
        
        sites_array = np.array(self.sites)
        self.ax.scatter(sites_array[:, 0], sites_array[:, 1], 
                       c=self.colors['sites'], s=80, zorder=5, 
                       edgecolors='white', linewidth=2, label='Sites')
    
    def plot_obstacles(self):
        """Plot obstacle points"""
        if not self.show_obstacles or not self.obstacles:
            return
        
        obstacles_array = np.array(self.obstacles)
        self.ax.scatter(obstacles_array[:, 0], obstacles_array[:, 1], 
                       c=self.colors['obstacles'], s=100, marker='X', 
                       zorder=5, edgecolors='white', linewidth=2, 
                       label='Obstacles')
    
    def update_plot(self):
        """Update the entire plot"""
        self.ax.clear()
        self.setup_plot()
        
        # Plot all components
        self.plot_voronoi()
        self.plot_convex_hull()
        self.plot_largest_empty_circle()
        self.plot_sites()
        self.plot_obstacles()
        
        # Add legend if there are any elements
        labels_present = []
        if self.sites:
            labels_present.append('Sites')
        if self.obstacles:
            labels_present.append('Obstacles')
        if self.convex_hull:
            labels_present.append('Convex Hull')
        if self.largest_empty_circle:
            labels_present.append('Largest Empty Circle')
        
        if labels_present:
            self.ax.legend(loc='upper right', framealpha=0.9)
        
        # Add statistics
        self.add_statistics()
        
        plt.draw()
    
    def add_statistics(self):
        """Add statistics text to the plot"""
        stats_text = []
        stats_text.append(f"Sites: {len(self.sites)}")
        stats_text.append(f"Obstacles: {len(self.obstacles)}")
        
        if self.voronoi:
            stats_text.append(f"Voronoi Edges: {len(self.voronoi.ridge_vertices)}")
        
        if self.convex_hull:
            stats_text.append(f"Hull Vertices: {len(self.convex_hull.vertices)}")
        
        if self.largest_empty_circle:
            radius = self.largest_empty_circle['radius']
            stats_text.append(f"LEC Radius: {radius:.3f}")
        
        # Add text box
        textstr = '\n'.join(stats_text)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        self.ax.text(0.02, 0.98, textstr, transform=self.ax.transAxes, fontsize=10,
                    verticalalignment='top', bbox=props)
    
    def on_click(self, event):
        """Handle mouse clicks to add/remove points"""
        if event.inaxes != self.ax:
            return
        
        # Check for modifier keys
        ctrl_pressed = event.key == 'control' if event.key else False
        
        # Add site point on left click
        if event.button == 1 and not ctrl_pressed:  # Left click
            self.add_site(event.xdata, event.ydata)
        # Add obstacle on right click
        elif event.button == 3:  # Right click
            self.add_obstacle(event.xdata, event.ydata)
        # Remove nearest point on middle click or Ctrl+Left click
        elif event.button == 2 or (event.button == 1 and ctrl_pressed):
            self.remove_nearest_point(event.xdata, event.ydata)
    
    def on_key_press(self, event):
        """Handle keyboard shortcuts"""
        key = event.key.lower() if event.key else ''
        
        if key == 'c':
            self.clear_all()
        elif key == 'r':
            self.generate_random_points()
        elif key == 'v':
            self.show_voronoi = not self.show_voronoi
            print(f"Voronoi diagram: {'ON' if self.show_voronoi else 'OFF'}")
            self.update_plot()
        elif key == 'h':
            self.show_convex_hull = not self.show_convex_hull
            print(f"Convex hull: {'ON' if self.show_convex_hull else 'OFF'}")
            self.update_plot()
        elif key == 'e':
            self.show_largest_circle = not self.show_largest_circle
            print(f"Largest empty circle: {'ON' if self.show_largest_circle else 'OFF'}")
            self.update_plot()
        elif key == 's':
            self.save_plot()
        elif key == 'help':
            self.setup_event_handlers()  # Print help again
    
    def save_plot(self, filename="voronoi_diagram.png", dpi=300):
        """Save the current plot"""
        self.fig.savefig(filename, dpi=dpi, bbox_inches='tight', 
                        facecolor='white', edgecolor='none')
        print(f"Plot saved as {filename}")
    
    def show(self):
        """Display the plot"""
        plt.show()

# Simple function to start the interactive plotter
def start_interactive_plotter():
    """Start the interactive Voronoi plotter with empty graph"""
    plotter = InteractiveVoronoiPlotter()
    plotter.show()
    return plotter

# Main execution
if __name__ == "__main__":
    # Start with empty interactive plotter
    plotter = start_interactive_plotter()