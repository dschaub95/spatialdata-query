import numpy as np

# Configure random seed for reproducibility
RNG = np.random.default_rng(2021)

def generate_points(n_points, bbox=(0, 100)):
    """Generate random points within the specified bounding box"""
    minx, maxx = bbox
    return RNG.random((int(n_points), 2)) * (maxx - minx) + minx

def generate_polygons(n_points, length=100):
    """Generate simple polygons - just squares for simplicity"""
    polygons = []
    for _ in range(int(n_points)):
        # Generate a random square
        x = RNG.random() * length
        y = RNG.random() * length
        size = RNG.random() * 10 + 1  # Size between 1 and 11
        
        # Square coordinates (5 points to close the polygon)
        polygon = [
            [x, y],
            [x + size, y],
            [x + size, y + size],
            [x, y + size],
            [x, y]  # Close the polygon
        ]
        polygons.append(polygon)
    
    return polygons