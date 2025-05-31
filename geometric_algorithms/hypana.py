import numpy as np
import pandas as pd
import scipy.interpolate
import scipy.ndimage as ndi
from scipy.optimize import minimize_scalar
# from skimage import measure
# import plotly.graph_objects as go
# import plotly.express as px
from typing import List, Dict, Tuple, Optional
import warnings
# import ripser
# import kmapper as km
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans



class ReactionHyperspaceAnalyzer:
    """
    Geometric analysis of reaction condition spaces for chemical optimization.

    Designed for experimental datasets where reaction conditions (temperature,
    concentrations, etc.) map to reaction yields. Handles both 3D and 4D
    hyperspace analysis with flexible column naming.

    Example usage:
        # For 3D analysis
        analyzer = ReactionHyperspaceAnalyzer(
            csv_path='hantzsch_data.csv',
            coordinate_columns=['temperature', 'conc_aldehyde', 'conc_amine'],
            yield_column='product_yield'
        )

        # For 4D analysis
        analyzer = ReactionHyperspaceAnalyzer(
            csv_path='ugi_data.csv',
            coordinate_columns=['conc_aldehyde', 'conc_amine', 'conc_isocyanide', 'conc_ptsa'],
            yield_column='yield_16e'
        )
    """

    def __init__(self, csv_path: str, coordinate_columns: List[str], yield_column: str):
        """
        Initialize analyzer with flexible column specification.

        Parameters:
        -----------
        csv_path : str
            Path to CSV file containing experimental data
        coordinate_columns : List[str]
            Names of columns containing reaction conditions (3 or 4 columns)
            Example: ['temperature', 'conc_A', 'conc_B'] for 3D
        yield_column : str
            Name of column containing reaction yields (0.0 to 1.0)
            Example: 'product_yield' or 'yield'
        """
        self.csv_path = csv_path
        self.coordinate_columns = coordinate_columns
        self.yield_column = yield_column
        self.dimensionality = len(coordinate_columns)

        # Validate inputs
        if self.dimensionality not in [3, 4]:
            raise ValueError(f"Only 3D and 4D analysis supported. Got {self.dimensionality} coordinates.")

        # Load and validate data
        self._load_and_validate_data()

        # Will be populated by interpolation
        self.regular_grid = None
        self.grid_coordinates = None
        self.yield_array = None
        self.interpolation_quality_metrics = {}

    def _load_and_validate_data(self):
        """Load CSV data and validate column specifications."""
        print(f"Loading data from {self.csv_path}...")

        try:
            self.raw_data = pd.read_csv(self.csv_path)
        except Exception as e:
            raise FileNotFoundError(f"Could not load CSV file: {e}")

        # Check that specified columns exist
        missing_cols = []
        for col in self.coordinate_columns + [self.yield_column]:
            if col not in self.raw_data.columns:
                missing_cols.append(col)

        if missing_cols:
            available_cols = list(self.raw_data.columns)
            raise ValueError(
                f"Missing columns: {missing_cols}\n"
                f"Available columns: {available_cols}"
            )

        # replace negative yields with zero
        self.raw_data[self.yield_column] = self.raw_data[self.yield_column].clip(lower=0.0)

        # Extract coordinate and yield data
        self.coordinates = self.raw_data[self.coordinate_columns].values
        self.yields = self.raw_data[self.yield_column].values

        # Basic data validation
        if np.any(np.isnan(self.coordinates)) or np.any(np.isnan(self.yields)):
            warnings.warn("Data contains NaN values. Consider cleaning data first.")

        if not (0 <= self.yields.min() and self.yields.max() <= 1.2):
            warnings.warn(
                f"Yields outside expected range [0,1]: min={self.yields.min():.3f}, max={self.yields.max():.3f}")

        print(f"Loaded {len(self.raw_data)} data points for {self.dimensionality}D analysis")
        print(f"Coordinate ranges:")
        for i, col in enumerate(self.coordinate_columns):
            coord_min, coord_max = self.coordinates[:, i].min(), self.coordinates[:, i].max()
            print(f"  {col}: [{coord_min:.3f}, {coord_max:.3f}]")
        print(f"Yield range: [{self.yields.min():.3f}, {self.yields.max():.3f}]")

    def interpolate_to_regular_grid(self, grid_resolution: int = 50,
                                    interpolation_method: str = 'linear'):
        """
        Resample irregular experimental data to regular N-dimensional grid.

        This is essential for applying geometric analysis methods that assume
        regular grid structure. Uses scipy's LinearNDInterpolator for robust
        interpolation of scattered data.

        Parameters:
        -----------
        grid_resolution : int, default=50
            Number of points along each coordinate axis
        interpolation_method : str, default='linear'
            Interpolation method ('linear', 'nearest', 'cubic' where supported)
        """
        print(f"Interpolating to regular {grid_resolution}^{self.dimensionality} grid...")

        # Create regular grid coordinates
        grid_axes = []
        self.coordinate_ranges = []

        for i, col_name in enumerate(self.coordinate_columns):
            coord_min = self.coordinates[:, i].min()
            coord_max = self.coordinates[:, i].max()

            # Add small padding to avoid edge effects
            padding = (coord_max - coord_min) * 0.02
            coord_min -= padding
            coord_max += padding

            axis = np.linspace(coord_min, coord_max, grid_resolution)
            grid_axes.append(axis)
            self.coordinate_ranges.append((coord_min, coord_max))

        # Create meshgrid for regular sampling points
        if self.dimensionality == 3:
            X, Y, Z = np.meshgrid(*grid_axes, indexing='ij')
            regular_points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
        elif self.dimensionality == 4:
            X, Y, Z, W = np.meshgrid(*grid_axes, indexing='ij')
            regular_points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel(), W.ravel()])

        # Perform interpolation
        from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

        # Try linear interpolation first
        try:
            interpolator = LinearNDInterpolator(self.coordinates, self.yields, fill_value=np.nan)
            interpolated_yields = interpolator(regular_points)

            # Check for too many NaN values (indicates sparse data)
            nan_fraction = np.isnan(interpolated_yields).mean()
            if nan_fraction > 0.3:
                print(f"Warning: {nan_fraction:.1%} of interpolated points are NaN (sparse data)")
                print("Filling NaN values with nearest neighbor interpolation...")

                # Fill NaN values with nearest neighbor
                nn_interpolator = NearestNDInterpolator(self.coordinates, self.yields)
                nan_mask = np.isnan(interpolated_yields)
                interpolated_yields[nan_mask] = nn_interpolator(regular_points[nan_mask])

        except Exception as e:
            print(f"Linear interpolation failed: {e}")
            print("Falling back to nearest neighbor interpolation...")
            interpolator = NearestNDInterpolator(self.coordinates, self.yields)
            interpolated_yields = interpolator(regular_points)

        # Reshape to grid
        grid_shape = tuple([grid_resolution] * self.dimensionality)
        self.yield_array = interpolated_yields.reshape(grid_shape)
        self.grid_coordinates = grid_axes
        self.grid_resolution = grid_resolution

        # Store interpolation quality metrics
        self._assess_interpolation_quality()

        print(f"Interpolation complete. Grid shape: {self.yield_array.shape}")
        print(f"Interpolated yield range: [{np.nanmin(self.yield_array):.3f}, {np.nanmax(self.yield_array):.3f}]")

    def save_regular_grid_to_file(self, filename: str):
        """
        Save interpolated regular grid to file for fast reloading.

        Saves the yield array and all necessary metadata to avoid
        re-interpolation of large 4D datasets.

        Parameters:
        -----------
        filename : str
            Output filename (will save as .npz format)
        """
        if self.yield_array is None:
            raise ValueError("No regular grid to save. Call interpolate_to_regular_grid() first.")

        # Ensure .npz extension
        if not filename.endswith('.npz'):
            filename += '.npz'

        print(f"Saving regular grid to {filename}...")

        # Save all necessary data to reconstruct the grid state
        np.savez_compressed(
            filename,
            yield_array=self.yield_array,
            grid_coordinates_0=self.grid_coordinates[0],
            grid_coordinates_1=self.grid_coordinates[1],
            grid_coordinates_2=self.grid_coordinates[2] if self.dimensionality >= 3 else np.array([]),
            grid_coordinates_3=self.grid_coordinates[3] if self.dimensionality >= 4 else np.array([]),
            coordinate_ranges=np.array(self.coordinate_ranges),
            grid_resolution=self.grid_resolution,
            dimensionality=self.dimensionality,
            coordinate_columns=np.array(self.coordinate_columns, dtype='U'),
            yield_column=np.array([self.yield_column], dtype='U'),
            interpolation_quality_metrics=np.array([str(self.interpolation_quality_metrics)], dtype='U')
        )

        print(f"Grid saved successfully. Shape: {self.yield_array.shape}")

    def load_regular_grid_from_file(self, filename: str):
        """
        Load previously saved regular grid from file.

        Reconstructs the interpolated grid state, bypassing the time-intensive
        interpolation step for large 4D datasets.

        Parameters:
        -----------
        filename : str
            Input filename (.npz format)
        """
        # Ensure .npz extension
        if not filename.endswith('.npz'):
            filename += '.npz'

        print(f"Loading regular grid from {filename}...")

        try:
            data = np.load(filename)
        except FileNotFoundError:
            raise FileNotFoundError(f"Grid file not found: {filename}")
        except Exception as e:
            raise ValueError(f"Could not load grid file: {e}")

        # Validate that the loaded grid matches current object configuration
        loaded_dimensionality = int(data['dimensionality'])
        loaded_coord_cols = data['coordinate_columns'].tolist()
        loaded_yield_col = data['yield_column'][0]

        if loaded_dimensionality != self.dimensionality:
            raise ValueError(f"Dimensionality mismatch: current={self.dimensionality}, loaded={loaded_dimensionality}")

        if loaded_coord_cols != self.coordinate_columns:
            raise ValueError(
                f"Coordinate columns mismatch: current={self.coordinate_columns}, loaded={loaded_coord_cols}")

        if loaded_yield_col != self.yield_column:
            raise ValueError(f"Yield column mismatch: current={self.yield_column}, loaded={loaded_yield_col}")

        # Reconstruct grid state
        self.yield_array = data['yield_array']
        self.grid_resolution = int(data['grid_resolution'])
        self.coordinate_ranges = data['coordinate_ranges'].tolist()

        # Reconstruct grid_coordinates list
        self.grid_coordinates = []
        for i in range(self.dimensionality):
            coord_data = data[f'grid_coordinates_{i}']
            if len(coord_data) > 0:  # Check for empty arrays (unused dimensions)
                self.grid_coordinates.append(coord_data)

        # Reconstruct interpolation quality metrics if available
        try:
            import ast
            metrics_str = data['interpolation_quality_metrics'][0]
            self.interpolation_quality_metrics = ast.literal_eval(metrics_str)
        except:
            self.interpolation_quality_metrics = {}

        print(f"Grid loaded successfully. Shape: {self.yield_array.shape}")
        print(f"Yield range: [{np.nanmin(self.yield_array):.3f}, {np.nanmax(self.yield_array):.3f}]")

    def _assess_interpolation_quality(self):
        """Assess quality of grid interpolation."""
        # Basic statistics
        self.interpolation_quality_metrics = {
            'original_points': len(self.yields),
            'grid_points': self.yield_array.size,
            'original_yield_range': (self.yields.min(), self.yields.max()),
            'interpolated_yield_range': (np.nanmin(self.yield_array), np.nanmax(self.yield_array)),
            'nan_fraction': np.isnan(self.yield_array).mean()
        }

        # Compute interpolation error at original data points
        if hasattr(self, 'yield_array'):
            from scipy.interpolate import RegularGridInterpolator

            # Create interpolator from regular grid back to original points
            grid_interpolator = RegularGridInterpolator(
                self.grid_coordinates, self.yield_array,
                method='linear', bounds_error=False, fill_value=np.nan
            )

            interpolated_at_original = grid_interpolator(self.coordinates)
            valid_mask = ~np.isnan(interpolated_at_original)

            if valid_mask.sum() > 10:  # Need sufficient points for meaningful statistics
                errors = np.abs(self.yields[valid_mask] - interpolated_at_original[valid_mask])
                self.interpolation_quality_metrics.update({
                    'mean_absolute_error': errors.mean(),
                    'max_absolute_error': errors.max(),
                    'rmse': np.sqrt(np.mean(errors ** 2)),
                    'validation_points': valid_mask.sum()
                })

                print(f"Interpolation quality (validated on {valid_mask.sum()} points):")
                print(f"  Mean absolute error: {errors.mean():.4f}")
                print(f"  RMS error: {np.sqrt(np.mean(errors ** 2)):.4f}")
                print(f"  Max error: {errors.max():.4f}")

    def find_local_maxima(self, min_distance: int = 3, threshold_abs: float = 0.1):
        """
        Identify local yield maxima with enforced minimum separation.

        Uses a two-step process:
        1. Find all candidate local maxima
        2. Greedily select peaks with minimum distance constraint

        Parameters:
        -----------
        min_distance : int, default=3
            Minimum distance between peaks (in grid points, Euclidean distance)
        threshold_abs : float, default=0.1
            Minimum absolute yield value for a peak

        Returns:
        --------
        Dict containing peak information
        """
        if self.yield_array is None:
            raise ValueError("Must call interpolate_to_regular_grid() first")

        print(f"Finding local maxima (min_distance={min_distance}, threshold={threshold_abs})...")

        # Step 1: Find all candidate local maxima using small neighborhood
        # Use small neighborhood (3^n) to find all local peaks first
        small_neighborhood = tuple([3] * self.dimensionality)
        candidate_maxima = (self.yield_array == ndi.maximum_filter(self.yield_array, size=small_neighborhood))

        # Apply threshold
        above_threshold = self.yield_array > threshold_abs
        candidate_mask = candidate_maxima & above_threshold

        # Get candidate indices and yields
        candidate_indices = np.argwhere(candidate_mask)

        if len(candidate_indices) == 0:
            print("No candidate maxima found. Try lowering threshold_abs.")
            return {
                'indices': np.array([]),
                'coordinates': np.array([]),
                'yields': np.array([]),
                'count': 0
            }

        candidate_yields = np.array([self.yield_array[tuple(idx)] for idx in candidate_indices])

        print(f"Found {len(candidate_indices)} candidate maxima before distance filtering...")

        # Step 2: Greedy selection with distance constraint
        # Sort candidates by yield (highest first)
        sort_order = np.argsort(candidate_yields)[::-1]
        sorted_indices = candidate_indices[sort_order]
        sorted_yields = candidate_yields[sort_order]

        # Greedily select peaks with minimum distance constraint
        selected_indices = []
        selected_yields = []

        for i, (idx, yield_val) in enumerate(zip(sorted_indices, sorted_yields)):
            # Check distance to all previously selected peaks
            too_close = False

            for selected_idx in selected_indices:
                # Calculate Euclidean distance in grid coordinates
                distance = np.sqrt(np.sum((idx - selected_idx) ** 2))

                if distance < min_distance:
                    too_close = True
                    break

            # If not too close to any existing peak, accept it
            if not too_close:
                selected_indices.append(idx)
                selected_yields.append(yield_val)

        if len(selected_indices) == 0:
            print("No peaks satisfy distance constraint. Try reducing min_distance.")
            return {
                'indices': np.array([]),
                'coordinates': np.array([]),
                'yields': np.array([]),
                'count': 0
            }

        # Convert to numpy arrays
        final_indices = np.array(selected_indices)
        final_yields = np.array(selected_yields)

        # Convert grid indices to parameter coordinates
        final_coordinates = np.zeros((len(final_indices), self.dimensionality))

        for i, idx in enumerate(final_indices):
            for dim in range(self.dimensionality):
                coord_axis = self.grid_coordinates[dim]
                final_coordinates[i, dim] = coord_axis[idx[dim]]

        print(f"Selected {len(final_indices)} peaks after distance filtering:")
        for i, (coord, yield_val) in enumerate(zip(final_coordinates, final_yields)):
            coord_str = ", ".join([f"{self.coordinate_columns[j]}={coord[j]:.3f}"
                                   for j in range(self.dimensionality)])
            print(f"  Peak {i + 1}: yield={yield_val:.3f} at ({coord_str})")

        # Verify minimum distances (for debugging)
        if len(final_indices) > 1:
            min_actual_distance = float('inf')
            for i in range(len(final_indices)):
                for j in range(i + 1, len(final_indices)):
                    dist = np.sqrt(np.sum((final_indices[i] - final_indices[j]) ** 2))
                    min_actual_distance = min(min_actual_distance, dist)
            print(f"  Minimum distance between selected peaks: {min_actual_distance:.1f} grid points")

        return {
            'indices': final_indices,
            'coordinates': final_coordinates,
            'yields': final_yields,
            'count': len(final_indices)
        }

    def compute_yield_gradients(self):
        """
        Compute gradient magnitude across the yield landscape.

        Uses finite differences to estimate gradients, providing insight into
        the 'steepness' of the yield surface. This directly addresses reviewer
        questions about smoothness and |D_ij| values.

        Returns:
        --------
        Dict containing gradient information
        """
        if self.yield_array is None:
            raise ValueError("Must call interpolate_to_regular_grid() first")

        print("Computing yield gradients...")

        # Compute gradients using numpy.gradient (central differences)
        gradients_tuple = np.gradient(self.yield_array)

        # Convert tuple to list for mutability
        gradients = list(gradients_tuple)

        # Convert to real coordinate units (not grid units)
        grid_spacings = []
        for i in range(self.dimensionality):
            coord_range = self.coordinate_ranges[i][1] - self.coordinate_ranges[i][0]
            grid_spacing = coord_range / (self.grid_resolution - 1)
            grid_spacings.append(grid_spacing)
            gradients[i] = gradients[i] / grid_spacing  # Now this works!

        # Compute gradient magnitude
        gradient_magnitude = np.sqrt(sum(g ** 2 for g in gradients))

        # Statistics (exclude NaN values)
        valid_mask = ~np.isnan(gradient_magnitude)
        if valid_mask.sum() > 0:
            mean_grad = np.mean(gradient_magnitude[valid_mask])
            max_grad = np.max(gradient_magnitude[valid_mask])

            print(f"Gradient statistics:")
            print(f"  Mean gradient magnitude: {mean_grad:.4f}")
            print(f"  Max gradient magnitude: {max_grad:.4f}")
            print(f"  Grid spacings: {[f'{s:.4f}' for s in grid_spacings]}")
        else:
            mean_grad = max_grad = np.nan
            print("Warning: All gradient values are NaN")

        return {
            'gradient_magnitude': gradient_magnitude,
            'gradient_components': gradients,
            'mean_gradient': mean_grad,
            'max_gradient': max_grad,
            'grid_spacings': grid_spacings
        }


    def validate_grid_quality(self):
        """
        Comprehensive validation of interpolated grid quality.

        Checks for common issues like sparse data coverage, interpolation
        artifacts, and provides recommendations for grid resolution.
        """
        if self.yield_array is None:
            raise ValueError("Must call interpolate_to_regular_grid() first")

        print("\n" + "=" * 50)
        print("GRID QUALITY VALIDATION")
        print("=" * 50)

        # 1. Data coverage
        nan_fraction = np.isnan(self.yield_array).mean()
        print(f"Data coverage: {(1 - nan_fraction) * 100:.1f}% of grid points have valid data")

        if nan_fraction > 0.5:
            print("⚠️  WARNING: More than 50% of grid points are NaN (very sparse data)")
        elif nan_fraction > 0.2:
            print("⚠️  CAUTION: More than 20% of grid points are NaN (somewhat sparse data)")
        else:
            print("✅ Good data coverage")

        # 2. Interpolation accuracy
        if 'mean_absolute_error' in self.interpolation_quality_metrics:
            mae = self.interpolation_quality_metrics['mean_absolute_error']
            yield_range = self.yields.max() - self.yields.min()
            relative_error = mae / yield_range * 100

            print(f"Interpolation accuracy: {relative_error:.1f}% relative error")
            if relative_error < 5:
                print("✅ Excellent interpolation accuracy")
            elif relative_error < 10:
                print("✅ Good interpolation accuracy")
            else:
                print("⚠️  WARNING: Poor interpolation accuracy - consider denser experimental sampling")

        # 3. Yield range preservation
        original_range = self.yields.max() - self.yields.min()
        interpolated_range = np.nanmax(self.yield_array) - np.nanmin(self.yield_array)
        range_preservation = interpolated_range / original_range * 100

        print(f"Yield range preservation: {range_preservation:.1f}%")
        if range_preservation > 95:
            print("✅ Excellent range preservation")
        elif range_preservation > 85:
            print("✅ Good range preservation")
        else:
            print("⚠️  CAUTION: Significant yield range compression during interpolation")

        # 4. Recommendations
        print("\nRecommendations:")
        if nan_fraction > 0.3:
            print("- Consider increasing experimental data density")
            print("- Try nearest neighbor interpolation if linear fails")

        if hasattr(self,
                   'interpolation_quality_metrics') and 'mean_absolute_error' in self.interpolation_quality_metrics:
            if self.interpolation_quality_metrics['mean_absolute_error'] > 0.05:
                print("- Consider increasing grid resolution")
                print("- Validate that experimental noise level is acceptable")

        print("=" * 50 + "\n")


    def analyze_maxima_basins(self, maxima_result: Dict, basin_threshold_fraction: float = 0.8):
        """
        Analyze local properties around each identified maximum.

        Computes curvature, basin volume, and robustness metrics for each peak.
        This provides quantitative characterization of the reaction landscape
        geometry around optimal conditions.

        Parameters:
        -----------
        maxima_result : Dict
            Output from find_local_maxima()
        basin_threshold_fraction : float, default=0.8
            Fraction of maximum yield to define basin boundary

        Returns:
        --------
        Dict with detailed analysis for each maximum
        """
        if self.yield_array is None:
            raise ValueError("Must call interpolate_to_regular_grid() first")

        if maxima_result['count'] == 0:
            print("No maxima to analyze")
            return {'basins': []}

        print(f"Analyzing {maxima_result['count']} maxima basins...")

        basin_analyses = []

        for i, (idx, coord, yield_val) in enumerate(zip(
                maxima_result['indices'],
                maxima_result['coordinates'],
                maxima_result['yields']
        )):
            print(f"Analyzing basin {i + 1}/{maxima_result['count']}...")

            basin_analysis = {
                'peak_index': i,
                'peak_coordinates': coord,
                'peak_yield': yield_val,
                'grid_indices': idx
            }

            # 1. Compute local Hessian matrix (curvature)
            hessian_info = self._compute_local_hessian(idx)
            basin_analysis.update(hessian_info)

            # 2. Estimate basin volume and robustness
            basin_info = self._analyze_basin_volume(idx, yield_val, basin_threshold_fraction)
            basin_analysis.update(basin_info)

            # 3. Compute gradient decay around maximum
            gradient_info = self._analyze_gradient_decay(idx, yield_val)
            basin_analysis.update(gradient_info)

            basin_analyses.append(basin_analysis)

        # 4. Comparative analysis between basins
        comparative_analysis = self._compare_basins(basin_analyses)

        return {
            'basins': basin_analyses,
            'comparative_analysis': comparative_analysis
        }


    def _compute_local_hessian(self, peak_idx: np.ndarray, neighborhood_size: int = 3):
        """Compute Hessian matrix around a peak for curvature analysis."""

        # Define neighborhood around peak
        slices = []
        for dim in range(self.dimensionality):
            start = max(0, peak_idx[dim] - neighborhood_size)
            end = min(self.yield_array.shape[dim], peak_idx[dim] + neighborhood_size + 1)
            slices.append(slice(start, end))

        # Extract local data
        local_yield = self.yield_array[tuple(slices)]

        # Compute second derivatives using finite differences
        hessian = np.zeros((self.dimensionality, self.dimensionality))

        # Get grid spacings in real coordinates
        grid_spacings = []
        for dim in range(self.dimensionality):
            coord_range = self.coordinate_ranges[dim][1] - self.coordinate_ranges[dim][0]
            spacing = coord_range / (self.grid_resolution - 1)
            grid_spacings.append(spacing)

        # Compute diagonal elements (pure second derivatives)
        for dim in range(self.dimensionality):
            if local_yield.shape[dim] >= 3:
                # Central difference for second derivative
                axis_slice = [slice(None)] * self.dimensionality

                # Get three points along this dimension
                center_idx = neighborhood_size
                if center_idx < local_yield.shape[dim] - 1 and center_idx > 0:
                    axis_slice[dim] = slice(center_idx - 1, center_idx + 2)
                    local_data = local_yield[tuple(axis_slice)]

                    # Extract values along the axis
                    if self.dimensionality == 3:
                        if dim == 0:
                            values = local_data[:, center_idx, center_idx]
                        elif dim == 1:
                            values = local_data[center_idx, :, center_idx]
                        else:
                            values = local_data[center_idx, center_idx, :]
                    elif self.dimensionality == 4:
                        remaining_idx = [center_idx] * 4
                        remaining_idx[dim] = slice(None)
                        # This is getting complex - use simpler approach
                        grad = np.gradient(local_yield, axis=dim)
                        grad2 = np.gradient(grad, axis=dim)
                        center_coord = tuple([center_idx if d != dim else slice(None) for d in range(self.dimensionality)])
                        if grad2[center_coord].size >= 3:
                            hessian[dim, dim] = grad2[center_coord][1] / (grid_spacings[dim] ** 2)
                        continue

                    if len(values) >= 3:
                        # Second derivative: (f(x+h) - 2*f(x) + f(x-h)) / h^2
                        h = grid_spacings[dim]
                        second_deriv = (values[2] - 2 * values[1] + values[0]) / (h ** 2)
                        hessian[dim, dim] = second_deriv

        # Compute eigenvalues for curvature characterization
        try:
            eigenvalues = np.linalg.eigvals(hessian)
            eigenvalues = eigenvalues[~np.isnan(eigenvalues)]

            if len(eigenvalues) > 0:
                max_curvature = np.max(np.abs(eigenvalues))
                min_curvature = np.min(eigenvalues)
                curvature_anisotropy = np.std(eigenvalues) if len(eigenvalues) > 1 else 0
            else:
                max_curvature = min_curvature = curvature_anisotropy = np.nan
        except:
            max_curvature = min_curvature = curvature_anisotropy = np.nan

        return {
            'hessian_matrix': hessian,
            'max_curvature': max_curvature,
            'min_curvature': min_curvature,
            'curvature_anisotropy': curvature_anisotropy,
            'hessian_eigenvalues': eigenvalues if 'eigenvalues' in locals() else np.array([])
        }


    def _analyze_basin_volume(self, peak_idx: np.ndarray, peak_yield: float, threshold_fraction: float):
        """Estimate basin volume and robustness around a maximum."""

        threshold_yield = peak_yield * threshold_fraction

        # Find connected component of high-yield region around this peak
        high_yield_mask = self.yield_array >= threshold_yield

        # Label connected components
        from scipy.ndimage import label
        labeled_array, num_features = label(high_yield_mask)

        # Find which component contains our peak
        peak_label = labeled_array[tuple(peak_idx)]

        if peak_label == 0:
            # Peak is not in high-yield region (shouldn't happen)
            basin_volume = 0
            basin_robustness = 0
        else:
            # Count voxels in the same component as our peak
            basin_mask = (labeled_array == peak_label)
            basin_volume = np.sum(basin_mask)

            # Compute robustness as average yield in basin in relation to peak yield
            basin_yields = self.yield_array[basin_mask]
            basin_robustness = np.mean(basin_yields) / peak_yield

        # Convert volume to real coordinate units
        voxel_volume = 1.0
        for dim in range(self.dimensionality):
            coord_range = self.coordinate_ranges[dim][1] - self.coordinate_ranges[dim][0]
            spacing = coord_range / (self.grid_resolution - 1)
            voxel_volume *= spacing

        basin_volume_real = basin_volume * voxel_volume

        return {
            'basin_volume_voxels': basin_volume,
            'basin_volume_real': basin_volume_real,
            'basin_robustness': basin_robustness,
            'threshold_used': threshold_yield
        }


    def _analyze_gradient_decay(self, peak_idx: np.ndarray, peak_yield: float, max_distance: int = 10):
        """Analyze how quickly yield drops off around the maximum."""

        distances = []
        yield_drops = []

        # Sample points at increasing distances from peak
        for distance in range(1, min(max_distance, self.grid_resolution // 4)):
            # Sample points on sphere/hypersphere around peak
            sampled_yields = []

            # Simple sampling: check points along each axis
            for dim in range(self.dimensionality):
                for direction in [-1, 1]:
                    test_idx = peak_idx.copy()
                    test_idx[dim] += direction * distance

                    # Check bounds
                    if 0 <= test_idx[dim] < self.yield_array.shape[dim]:
                        test_yield = self.yield_array[tuple(test_idx)]
                        if not np.isnan(test_yield):
                            sampled_yields.append(test_yield)

            if sampled_yields:
                avg_yield_at_distance = np.mean(sampled_yields)
                yield_drop = peak_yield - avg_yield_at_distance

                distances.append(distance)
                yield_drops.append(yield_drop)

        # Fit exponential decay if we have enough points
        decay_rate = np.nan
        if len(distances) >= 3:
            try:
                from scipy.optimize import curve_fit

                def exponential_decay(x, a, b):
                    return a * (1 - np.exp(-b * x))

                popt, _ = curve_fit(exponential_decay, distances, yield_drops,
                                    bounds=([0, 0], [1, 10]), maxfev=1000)
                decay_rate = popt[1]
            except:
                pass

        return {
            'gradient_decay_distances': np.array(distances),
            'gradient_decay_yield_drops': np.array(yield_drops),
            'exponential_decay_rate': decay_rate
        }


    def _compare_basins(self, basin_analyses: List[Dict]):
        """Compare properties between different basins."""

        if len(basin_analyses) <= 1:
            return {'note': 'Need at least 2 basins for comparison'}

        # Extract key metrics
        peak_yields = [b['peak_yield'] for b in basin_analyses]
        basin_volumes = [b['basin_volume_real'] for b in basin_analyses]
        max_curvatures = [b['max_curvature'] for b in basin_analyses if not np.isnan(b['max_curvature'])]
        decay_rates = [b['exponential_decay_rate'] for b in basin_analyses if not np.isnan(b['exponential_decay_rate'])]

        comparison = {
            'yield_ratio': max(peak_yields) / min(peak_yields) if min(peak_yields) > 0 else np.inf,
            'volume_ratio': max(basin_volumes) / min(basin_volumes) if min(basin_volumes) > 0 else np.inf,
            'curvature_variation': np.std(max_curvatures) if len(max_curvatures) > 1 else 0,
            'decay_rate_variation': np.std(decay_rates) if len(decay_rates) > 1 else 0
        }

        return comparison


    # Tier 2 implementations
    def compute_persistent_homology(self, superlevel: bool = True, max_dimension: int = 1):
        """
        Compute persistent homology to identify robust features in yield landscape.

        Uses ripser to find topological features like robust maxima (0-dimensional)
        and potential "bridges" between high-yield regions (1-dimensional).

        Parameters:
        -----------
        superlevel : bool, default=True
            If True, analyze superlevel sets (high yield regions)
            If False, analyze sublevel sets (low yield regions)
        max_dimension : int, default=1
            Maximum homological dimension to compute

        Returns:
        --------
        Dict containing persistence diagrams and analysis
        """
        if self.yield_array is None:
            raise ValueError("Must call interpolate_to_regular_grid() first")

        print("Computing persistent homology...")

        try:
            import ripser
        except ImportError:
            raise ImportError("ripser package required. Install with: pip install ripser")

        # Prepare data for ripser (needs 1D distance matrix or point cloud)
        # For cubical/grid data, we'll use the ripser's built-in cubical complex functionality

        # Flatten yield array and remove NaN values
        yield_flat = self.yield_array.ravel()
        valid_mask = ~np.isnan(yield_flat)

        if not valid_mask.any():
            raise ValueError("No valid yield data for persistence computation")

        # For superlevel sets, we need to negate the function
        if superlevel:
            yield_for_ripser = -yield_flat[valid_mask]
        else:
            yield_for_ripser = yield_flat[valid_mask]

        # Compute persistence using ripser
        # Note: ripser expects distance matrices or point clouds, so we'll use a workaround
        # Create point cloud in yield space for 0D persistence

        # Simple approach: treat each yield value as a 1D point
        points_1d = yield_for_ripser.reshape(-1, 1)

        try:
            result = ripser.ripser(points_1d, maxdim=max_dimension, distance_matrix=False)
            diagrams = result['dgms']

            print(f"Computed {len(diagrams)} persistence diagrams")
            for i, dgm in enumerate(diagrams):
                if len(dgm) > 0:
                    print(f"  Dimension {i}: {len(dgm)} features")
                    if i == 0:  # 0-dimensional features (connected components)
                        # Find long-lived features
                        persistences = dgm[:, 1] - dgm[:, 0]
                        long_lived = persistences > np.percentile(persistences, 90)
                        print(f"    {np.sum(long_lived)} robust features (top 10% persistence)")

            return {
                'diagrams': diagrams,
                'superlevel': superlevel,
                'long_lived_features_0d': self._extract_long_lived_features(
                    diagrams[0] if len(diagrams) > 0 else np.array([]))
            }

        except Exception as e:
            print(f"Ripser computation failed: {e}")
            return {'diagrams': [], 'error': str(e)}


    def _extract_long_lived_features(self, diagram_0d: np.ndarray, persistence_threshold: float = None):
        """Extract long-lived 0-dimensional features (robust maxima)."""

        if len(diagram_0d) == 0:
            return []

        # Compute persistence
        persistences = diagram_0d[:, 1] - diagram_0d[:, 0]

        # Auto-threshold if not provided
        if persistence_threshold is None:
            persistence_threshold = np.percentile(persistences, 90)

        # Find long-lived features
        long_lived_mask = persistences > persistence_threshold
        long_lived_features = diagram_0d[long_lived_mask]

        return {
            'features': long_lived_features,
            'count': len(long_lived_features),
            'persistence_threshold': persistence_threshold,
            'mean_persistence': np.mean(persistences[long_lived_mask]) if long_lived_mask.any() else 0
        }

    def create_mapper_graph(self, lens_function: str = 'yield', n_intervals: int = 15,
                            overlap: float = 0.3, clustering_algorithm: str = 'DBSCAN',
                            dbscan_eps: float = None, dbscan_min_samples: int = None):
        """
        Create Mapper graph for high-dimensional visualization.

        Particularly useful for 4D data where direct visualization is impossible.
        Creates a graph representation showing connectivity of high-yield regions.

        Parameters:
        -----------
        lens_function : str, default='yield'
            Function to use as lens ('yield', 'first_coordinate', etc.)
        n_intervals : int, default=15
            Number of intervals for cover
        overlap : float, default=0.3
            Overlap fraction between intervals
        clustering_algorithm : str, default='DBSCAN'
            Clustering method ('DBSCAN', 'KMeans')
        dbscan_eps : float, optional
            DBSCAN epsilon parameter. Auto-tuned if None.
        dbscan_min_samples : int, optional
            DBSCAN min_samples parameter. Auto-tuned if None.

        Returns:
        --------
        Dict containing mapper graph and metadata
        """
        if self.yield_array is None:
            raise ValueError("Must call interpolate_to_regular_grid() first")

        print(f"Creating Mapper graph with {lens_function} lens...")

        try:
            import kmapper as km
            from sklearn.cluster import DBSCAN
            from sklearn.cluster import KMeans
        except ImportError:
            raise ImportError("kmapper and scikit-learn required. Install with: pip install kmapper scikit-learn")

        # Prepare data: flatten grid to point cloud
        valid_mask = ~np.isnan(self.yield_array)
        valid_indices = np.argwhere(valid_mask)
        n_valid = len(valid_indices)

        if n_valid == 0:
            raise ValueError("No valid data points for Mapper")

        # Convert grid indices to real coordinates
        point_cloud = np.zeros((n_valid, self.dimensionality + 1))  # coords + yield

        for i, grid_idx in enumerate(valid_indices):
            # Convert grid index to real coordinates
            for dim in range(self.dimensionality):
                coord_value = self.grid_coordinates[dim][grid_idx[dim]]
                point_cloud[i, dim] = coord_value

            # Add yield value
            point_cloud[i, -1] = self.yield_array[tuple(grid_idx)]

        # Debug data properties
        print(f"Point cloud shape: {point_cloud.shape}")
        print(f"Yield range in point cloud: [{np.min(point_cloud[:, -1]):.4f}, {np.max(point_cloud[:, -1]):.4f}]")
        print(
            f"Coordinate ranges: {[(np.min(point_cloud[:, i]), np.max(point_cloud[:, i])) for i in range(self.dimensionality)]}")

        # CRITICAL: More aggressive subsampling for meaningful clustering
        target_points = min(1000, n_valid)  # Much smaller for better clustering
        if n_valid > target_points:
            subsample_idx = np.random.choice(n_valid, target_points, replace=False)
            point_cloud = point_cloud[subsample_idx]
            print(f"Subsampled to {len(point_cloud)} points for Mapper")

        # Set up lens function
        if lens_function == 'yield':
            lens = point_cloud[:, -1].reshape(-1, 1)
        elif lens_function == 'first_coordinate':
            lens = point_cloud[:, 0].reshape(-1, 1)
        else:
            raise ValueError(f"Unknown lens function: {lens_function}")

        print(f"Lens range: [{np.min(lens):.4f}, {np.max(lens):.4f}]")

        # CRITICAL FIX: Use only coordinate data for clustering (not yield!)
        clustering_data = point_cloud[:, :-1]  # Exclude yield column

        # FIXED: Calculate coord_ranges before any conditional logic
        coord_ranges = np.max(clustering_data, axis=0) - np.min(clustering_data, axis=0)
        print(f"Coordinate ranges for clustering: {coord_ranges}")

        # Set up clustering with proper scaling
        if clustering_algorithm == 'DBSCAN':
            if dbscan_eps is None:
                # Scale eps to coordinate ranges
                dbscan_eps = np.mean(coord_ranges) * 0.02  # 2% of average coordinate range

            if dbscan_min_samples is None:
                dbscan_min_samples = max(3, len(point_cloud) // 100)  # At least 3, scale with data

            print(f"DBSCAN parameters: eps={dbscan_eps:.6f}, min_samples={dbscan_min_samples}")
            clusterer = DBSCAN(eps=dbscan_eps, min_samples=dbscan_min_samples)

            # Test clustering on small sample
            test_clustering = DBSCAN(eps=dbscan_eps, min_samples=dbscan_min_samples)
            test_labels = test_clustering.fit_predict(clustering_data[:100])
            n_clusters_test = len(set(test_labels)) - (1 if -1 in test_labels else 0)
            print(f"Test clustering on 100 points: {n_clusters_test} clusters")

            if n_clusters_test == 0:
                print("WARNING: Test clustering found 0 clusters. Adjusting parameters...")
                # Try progressively larger eps values
                for multiplier in [0.05, 0.1, 0.15, 0.2]:
                    dbscan_eps = np.mean(coord_ranges) * multiplier
                    dbscan_min_samples = max(2, dbscan_min_samples // 2)

                    test_clustering = DBSCAN(eps=dbscan_eps, min_samples=dbscan_min_samples)
                    test_labels = test_clustering.fit_predict(clustering_data[:100])
                    n_clusters_test = len(set(test_labels)) - (1 if -1 in test_labels else 0)

                    print(
                        f"  Trying eps={dbscan_eps:.6f}, min_samples={dbscan_min_samples}: {n_clusters_test} clusters")

                    if n_clusters_test > 0:
                        break

                if n_clusters_test == 0:
                    print("ERROR: Could not find working DBSCAN parameters. Switching to KMeans.")
                    clustering_algorithm = 'KMeans'
                else:
                    clusterer = DBSCAN(eps=dbscan_eps, min_samples=dbscan_min_samples)
                    print(f"Final DBSCAN: eps={dbscan_eps:.6f}, min_samples={dbscan_min_samples}")

        if clustering_algorithm == 'KMeans':
            n_clusters = min(8, max(2, len(point_cloud) // 100))  # Reasonable number of clusters
            clusterer = KMeans(n_clusters=n_clusters, random_state=42)
            print(f"Using KMeans with {n_clusters} clusters")

        # Create Mapper with fewer intervals for cleaner graph
        mapper = km.KeplerMapper(verbose=1)

        try:
            graph = mapper.map(lens=lens,
                               X=clustering_data,
                               clusterer=clusterer,
                               cover=km.Cover(n_cubes=max(6, n_intervals // 3), perc_overlap=overlap))
        except Exception as e:
            return {'error': f'Mapper creation failed: {str(e)}'}

        # Analyze graph properties
        n_nodes = len(graph['nodes'])
        n_edges = len(graph['links'])

        print(f"Mapper graph: {n_nodes} nodes, {n_edges} edges")

        # Target: 10-50 nodes for meaningful visualization
        if n_nodes > 100:
            print(f"⚠️  Warning: {n_nodes} nodes is too many for clear visualization.")
            print("   Consider: increasing dbscan_eps, increasing dbscan_min_samples, fewer n_intervals")
        elif n_nodes < 5:
            print(f"⚠️  Warning: {n_nodes} nodes may be too few. Consider decreasing dbscan_eps")
        else:
            print(f"✅ Good graph complexity: {n_nodes} nodes")

        # Find high-yield clusters
        high_yield_nodes = []
        for node_id, point_indices in graph['nodes'].items():
            node_yields = point_cloud[point_indices, -1]
            avg_yield = np.mean(node_yields)
            max_yield = np.max(node_yields)

            if avg_yield > 0.1:  # Threshold for your data
                high_yield_nodes.append({
                    'node_id': node_id,
                    'avg_yield': avg_yield,
                    'max_yield': max_yield,
                    'size': len(point_indices)
                })

        print(f"Found {len(high_yield_nodes)} high-yield nodes")

        return {
            'graph': graph,
            'point_cloud': point_cloud,
            'lens': lens,
            'high_yield_nodes': high_yield_nodes,
            'n_nodes': n_nodes,
            'n_edges': n_edges,
            'lens_function': lens_function
        }


    def plot_mapper_graph_summary(self, mapper_result: Dict, save_path: str = None):
        """
        Create summary visualization of Mapper graph.

        Shows graph structure with nodes colored by average yield.
        Handles both simple and complex graphs appropriately.
        """
        if 'error' in mapper_result:
            print(f"Cannot plot Mapper graph: {mapper_result['error']}")
            return

        try:
            import matplotlib.pyplot as plt
            import networkx as nx
        except ImportError:
            raise ImportError("matplotlib and networkx required for plotting")

        graph = mapper_result['graph']
        point_cloud = mapper_result['point_cloud']

        # If too many nodes, use simplified visualization
        if len(graph['nodes']) > 100:
            print(f"Graph too complex ({len(graph['nodes'])} nodes). Creating simplified visualization...")
            return self._plot_simplified_mapper(mapper_result, save_path)

        print(f"Debug: Graph has {len(graph['nodes'])} nodes")
        print(f"Debug: Node IDs: {list(graph['nodes'].keys())[:5]}...")  # Show first 5

        # Create NetworkX graph
        G = nx.Graph()

        # Add nodes with yield information
        node_yields = {}
        node_sizes = {}

        for node_id, point_indices in graph['nodes'].items():
            # Ensure point_indices is valid
            if len(point_indices) > 0:
                node_yield = np.mean(point_cloud[point_indices, -1])
                # Handle NaN values
                if np.isnan(node_yield):
                    node_yield = 0.0
                node_yields[node_id] = node_yield
                node_sizes[node_id] = len(point_indices)
                G.add_node(node_id)

        print(f"Debug: Added {len(G.nodes())} nodes to NetworkX graph")

        # Add edges - handle different kmapper edge formats
        edges_added = 0
        for edge in graph['links']:
            try:
                if isinstance(edge, (list, tuple)) and len(edge) >= 2:
                    node1, node2 = edge[0], edge[1]
                    # Only add edge if both nodes exist
                    if node1 in G.nodes() and node2 in G.nodes():
                        G.add_edge(node1, node2)
                        edges_added += 1
                elif isinstance(edge, dict):
                    # Handle dictionary format: {'source': ..., 'target': ...}
                    if 'source' in edge and 'target' in edge:
                        node1, node2 = edge['source'], edge['target']
                        if node1 in G.nodes() and node2 in G.nodes():
                            G.add_edge(node1, node2)
                            edges_added += 1
            except Exception as e:
                print(f"Warning: Could not process edge {edge}: {e}")
                continue

        print(f"Debug: Added {edges_added} edges to NetworkX graph")

        if len(G.nodes()) == 0:
            print("Error: No valid nodes in graph for visualization")
            return None

        # Create visualization
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        # Plot 1: Graph colored by yield
        try:
            pos = nx.spring_layout(G, seed=42, k=1, iterations=50)
        except:
            # Fallback to random layout if spring layout fails
            pos = nx.random_layout(G, seed=42)

        # Extract data for plotting - ensure proper ordering
        node_list = list(G.nodes())
        yields = [node_yields[node] for node in node_list]
        sizes = [node_sizes[node] * 50 for node in node_list]  # Scale for visibility

        # Debug the data
        print(f"Debug: yields range: {min(yields):.3f} to {max(yields):.3f}")
        print(f"Debug: sizes range: {min(sizes)} to {max(sizes)}")

        # Check for problematic values
        yields = np.array(yields)
        sizes = np.array(sizes)

        # Replace any remaining NaN or inf values
        yields = np.nan_to_num(yields, nan=0.0, posinf=1.0, neginf=0.0)
        sizes = np.nan_to_num(sizes, nan=50.0, posinf=500.0, neginf=50.0)

        # Extract positions
        x_pos = [pos[node][0] for node in node_list]
        y_pos = [pos[node][1] for node in node_list]

        # Create scatter plot with proper error handling
        try:
            scatter = ax1.scatter(x_pos, y_pos, c=yields, s=sizes,
                                  cmap='viridis', alpha=0.7, vmin=yields.min(), vmax=yields.max())
        except Exception as e:
            print(f"Scatter plot failed: {e}")
            # Fallback: simple plot without colors
            scatter = ax1.scatter(x_pos, y_pos, s=sizes, alpha=0.7, c='blue')

        # Draw edges
        edges_drawn = 0
        for edge in G.edges():
            try:
                if edge[0] in pos and edge[1] in pos:
                    x_vals = [pos[edge[0]][0], pos[edge[1]][0]]
                    y_vals = [pos[edge[0]][1], pos[edge[1]][1]]
                    ax1.plot(x_vals, y_vals, 'k-', alpha=0.3, linewidth=0.5)
                    edges_drawn += 1
            except Exception as e:
                print(f"Edge drawing failed for edge {edge}: {e}")
                continue

        print(f"Debug: Drew {edges_drawn} edges")

        ax1.set_title(f'Mapper Graph ({len(G.nodes())} nodes, {len(G.edges())} edges)')
        ax1.set_xlabel('Spring Layout X')
        ax1.set_ylabel('Spring Layout Y')

        # Add colorbar if scatter plot worked
        try:
            if 'scatter' in locals() and hasattr(scatter, 'get_cmap'):
                cbar1 = plt.colorbar(scatter, ax=ax1)
                cbar1.set_label('Average Yield')
        except:
            pass

        # Plot 2: Yield distribution across nodes
        try:
            ax2.hist(yields, bins=min(10, len(yields)), alpha=0.7, edgecolor='black')
            ax2.set_xlabel('Average Yield')
            ax2.set_ylabel('Number of Nodes')
            ax2.set_title('Yield Distribution Across Graph Nodes')

            if len(yields) > 0:
                mean_yield = np.mean(yields)
                ax2.axvline(mean_yield, color='red', linestyle='--',
                            label=f'Mean: {mean_yield:.3f}')
                ax2.legend()
        except Exception as e:
            print(f"Histogram plot failed: {e}")
            ax2.text(0.5, 0.5, 'Histogram failed', transform=ax2.transAxes, ha='center')

        plt.tight_layout()

        if save_path:
            try:
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"Mapper graph saved to {save_path}")
            except Exception as e:
                print(f"Saving failed: {e}")

        try:
            plt.show()
        except:
            print("Display failed, but figure created")

        return fig


    def _plot_simplified_mapper(self, mapper_result: Dict, save_path: str = None):
        """
        Simplified visualization for very complex graphs.

        When the Mapper graph has too many nodes (>100), this creates
        aggregate visualizations focusing on high-yield regions and
        overall data structure rather than individual node connectivity.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib required for plotting")

        graph = mapper_result['graph']
        point_cloud = mapper_result['point_cloud']
        high_yield_nodes = mapper_result['high_yield_nodes']

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

        # Plot 1: High-yield node analysis
        if high_yield_nodes:
            yields = [node['avg_yield'] for node in high_yield_nodes]
            sizes = [node['size'] for node in high_yield_nodes]

            scatter = ax1.scatter(range(len(yields)), yields, s=sizes, alpha=0.7, c=yields, cmap='viridis')
            ax1.set_xlabel('High-Yield Node Index')
            ax1.set_ylabel('Average Yield')
            ax1.set_title(f'High-Yield Nodes ({len(high_yield_nodes)} total)')

            try:
                plt.colorbar(scatter, ax=ax1, label='Average Yield')
            except:
                pass

            # Add text annotations for top nodes
            for i, (y, s) in enumerate(zip(yields[:5], sizes[:5])):  # Top 5 nodes
                ax1.annotate(f'{y:.3f}', (i, y), xytext=(5, 5), textcoords='offset points', fontsize=8)
        else:
            ax1.text(0.5, 0.5, 'No high-yield nodes found', transform=ax1.transAxes, ha='center')
            ax1.set_title('High-Yield Nodes (None Found)')

        # Plot 2: Overall yield distribution
        all_yields = point_cloud[:, -1]
        ax2.hist(all_yields, bins=30, alpha=0.7, edgecolor='black')
        ax2.set_xlabel('Yield')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Overall Yield Distribution')
        ax2.axvline(np.mean(all_yields), color='red', linestyle='--',
                    label=f'Mean: {np.mean(all_yields):.3f}')
        ax2.axvline(np.median(all_yields), color='orange', linestyle='--',
                    label=f'Median: {np.median(all_yields):.3f}')
        ax2.legend()

        # Plot 3: Node size distribution
        all_node_sizes = [len(point_indices) for point_indices in graph['nodes'].values()]
        ax3.hist(all_node_sizes, bins=20, alpha=0.7, edgecolor='black')
        ax3.set_xlabel('Points per Node')
        ax3.set_ylabel('Number of Nodes')
        ax3.set_title('Node Size Distribution')
        ax3.axvline(np.mean(all_node_sizes), color='red', linestyle='--',
                    label=f'Mean: {np.mean(all_node_sizes):.1f}')
        ax3.legend()

        # Plot 4: Graph complexity metrics
        graph_metrics = {
            'Total Nodes': len(graph['nodes']),
            'Total Edges': len(graph['links']),
            'High-Yield Nodes': len(high_yield_nodes),
            'Avg Points/Node': np.mean(all_node_sizes),
            'Edge Density': len(graph['links']) / max(1, len(graph['nodes']) * (len(graph['nodes']) - 1) / 2)
        }

        metrics_text = '\n'.join([f'{k}: {v:.3f}' if isinstance(v, float) else f'{k}: {v}'
                                  for k, v in graph_metrics.items()])

        ax4.text(0.1, 0.9, 'Graph Complexity Metrics:', transform=ax4.transAxes,
                 fontsize=12, fontweight='bold', verticalalignment='top')
        ax4.text(0.1, 0.7, metrics_text, transform=ax4.transAxes,
                 fontsize=10, verticalalignment='top', fontfamily='monospace')

        # Add interpretation
        interpretation = []
        if len(high_yield_nodes) >= 2:
            interpretation.append("✅ Multiple high-yield regions detected")
        elif len(high_yield_nodes) == 1:
            interpretation.append("⚠️  Single high-yield region found")
        else:
            interpretation.append("❌ No high-yield regions found")

        if len(graph['nodes']) > 200:
            interpretation.append("⚠️  Very complex graph - consider\n   increasing clustering parameters")
        elif len(graph['nodes']) < 10:
            interpretation.append("⚠️  Very simple graph - consider\n   decreasing clustering parameters")
        else:
            interpretation.append("✅ Reasonable graph complexity")

        ax4.text(0.1, 0.4, 'Interpretation:', transform=ax4.transAxes,
                 fontsize=12, fontweight='bold', verticalalignment='top')
        ax4.text(0.1, 0.3, '\n'.join(interpretation), transform=ax4.transAxes,
                 fontsize=10, verticalalignment='top')

        ax4.set_xlim(0, 1)
        ax4.set_ylim(0, 1)
        ax4.axis('off')

        plt.tight_layout()

        if save_path:
            try:
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"Simplified mapper plot saved to {save_path}")
            except Exception as e:
                print(f"Saving failed: {e}")

        try:
            plt.show()
        except:
            print("Display failed, but figure created")

        return fig


# Usage example for testing
if __name__ == "__main__":
    pass
    # # ### For 3D analysis
    # analyzer = ReactionHyperspaceAnalyzer(
    #     csv_path='simpleSN1_resampled_yields_hyvu.csv',
    #     coordinate_columns=['[Alcohol](mM)','[HBr](mM)','Temperature(°C)'],
    #     yield_column='yield'
    # )
    #
    # # Interpolate to regular grid
    # analyzer.interpolate_to_regular_grid(grid_resolution=5)
    #
    # analyzer.save_regular_grid_to_file('3d_grid.npz')
    #
    # analyzer.load_regular_grid_from_file('3d_grid.npz')
    #
    # # # Validate quality
    # # analyzer.validate_grid_quality()
    #
    # # Find local maxima
    # maxima = analyzer.find_local_maxima(min_distance=2, threshold_abs=0.9)
    #
    # # Compute gradients
    # gradients = analyzer.compute_yield_gradients()

    # For 4D analysis
    analyzer = ReactionHyperspaceAnalyzer(
        csv_path='postprocessed_yields_decimated.csv',
        coordinate_columns=['ic001','am001','ald001','ptsa'],
        yield_column='yield'
    )

    # Interpolate to regular grid
    analyzer.interpolate_to_regular_grid(grid_resolution=10)

    analyzer.save_regular_grid_to_file('4d_grid.npz')

    analyzer.load_regular_grid_from_file('4d_grid.npz')

    # Validate quality
    analyzer.validate_grid_quality()

    # Find local maxima
    maxima = analyzer.find_local_maxima(min_distance=5, threshold_abs=0.11)

    # # Compute gradients
    # gradients = analyzer.compute_yield_gradients()

    # Test the new functionality:

    print("\n" + "=" * 60)
    print("TESTING NEW GEOMETRIC ANALYSIS FUNCTIONALITY")
    print("=" * 60)

    # 1. Analyze basins around the identified maxima
    print("\n1. BASIN ANALYSIS")
    print("-" * 30)
    basin_analysis = analyzer.analyze_maxima_basins(maxima, basin_threshold_fraction=0.8)

    # Print summary of basin analysis
    print(f"\nBasin Analysis Summary:")
    print(f"Number of basins analyzed: {len(basin_analysis['basins'])}")

    for i, basin in enumerate(basin_analysis['basins']):
        print(f"\nBasin {i + 1}:")
        print(f"  Peak yield: {basin['peak_yield']:.4f}")
        print(f"  Peak coordinates: {dict(zip(analyzer.coordinate_columns, basin['peak_coordinates']))}")
        print(f"  Basin volume (real units): {basin['basin_volume_real']}")
        print(f"  Basin robustness: {basin['basin_robustness']:.3f}")
        print(f"  Max curvature: {basin['max_curvature']:.4f}")

        if not np.isnan(basin['exponential_decay_rate']):
            print(f"  Gradient decay rate: {basin['exponential_decay_rate']:.4f}")

    # Print comparative analysis
    if 'comparative_analysis' in basin_analysis:
        comp = basin_analysis['comparative_analysis']
        print(f"\nComparative Analysis:")
        if 'yield_ratio' in comp:
            print(f"  Yield ratio (max/min): {comp['yield_ratio']:.3f}")
        if 'volume_ratio' in comp:
            print(f"  Volume ratio (max/min): {comp['volume_ratio']:.3f}")

    # 2. Persistent homology analysis
    print("\n\n2. PERSISTENT HOMOLOGY ANALYSIS")
    print("-" * 40)
    ph_result = analyzer.compute_persistent_homology(superlevel=True, max_dimension=1)

    if 'error' not in ph_result:
        print("Persistent homology computed successfully!")

        # Analyze 0-dimensional features (should show your two maxima)
        if 'long_lived_features_0d' in ph_result and ph_result['long_lived_features_0d']:
            features_0d = ph_result['long_lived_features_0d']
            print(f"Long-lived 0D features (robust maxima): {features_0d['count']}")
            print(f"Mean persistence: {features_0d['mean_persistence']:.4f}")

            # This should show 2 features for your two-maxima dataset!
            if features_0d['count'] == 2:
                print("✅ SUCCESS: Found exactly 2 robust maxima (matches expected two disjoined maxima!)")
            elif features_0d['count'] == 1:
                print("⚠️  Found only 1 robust maximum - may need to adjust persistence threshold")
            else:
                print(f"⚠️  Found {features_0d['count']} robust maxima - unexpected for this dataset")

        # Print dimensions of persistence diagrams
        for i, dgm in enumerate(ph_result['diagrams']):
            print(f"Dimension {i} persistence: {len(dgm)} features")
    else:
        print(f"Persistent homology failed: {ph_result['error']}")

    # 3. Mapper graph analysis (great for 4D visualization)
    print("\n\n3. MAPPER GRAPH ANALYSIS (4D VISUALIZATION)")
    print("-" * 50)
    mapper_result = analyzer.create_mapper_graph(
        lens_function='yield',
        n_intervals=8,  # Fewer intervals
        overlap=0.4,  # Reasonable overlap
        clustering_algorithm='DBSCAN',
        dbscan_eps=0.01,  # Larger eps for bigger clusters
        dbscan_min_samples=10  # Much higher min_samples
    )

    if 'error' not in mapper_result:
        print("Mapper graph created successfully!")
        print(f"Graph structure: {mapper_result['n_nodes']} nodes, {mapper_result['n_edges']} edges")
        print(f"High-yield nodes: {len(mapper_result['high_yield_nodes'])}")

        # Analyze high-yield nodes (should cluster into 2 groups for your data)
        if mapper_result['high_yield_nodes']:
            high_yields = [node['avg_yield'] for node in mapper_result['high_yield_nodes']]
            print(f"High-yield node yields: {[f'{y:.3f}' for y in sorted(high_yields, reverse=True)]}")

            # Check if we can see separation (clustering of high-yield nodes)
            if len(high_yields) >= 2:
                yield_gap = max(high_yields) - min(high_yields)
                print(f"Yield gap between high-yield regions: {yield_gap:.3f}")

        # Create visualization
        print("\nCreating Mapper graph visualization...")
        try:
            analyzer.plot_mapper_graph_summary(mapper_result, save_path='mapper_4d_analysis.png')
        except Exception as e:
            print(f"Visualization failed: {e}")

    else:
        print(f"Mapper graph creation failed: {mapper_result['error']}")

    # 4. Overall gradient analysis (addresses reviewer smoothness questions)
    print("\n\n4. OVERALL GRADIENT ANALYSIS")
    print("-" * 35)
    gradient_analysis = analyzer.compute_yield_gradients()

    print(f"Mean gradient magnitude: {gradient_analysis['mean_gradient']:.6f}")
    print(f"Max gradient magnitude: {gradient_analysis['max_gradient']:.6f}")
    print(f"This addresses reviewer questions about |D_ij| ≤ 1 smoothness")

    # Check if gradients are indeed small (supporting your smoothness claims)
    if gradient_analysis['max_gradient'] <= 1.0:
        print("✅ Confirms smooth landscape (max gradient ≤ 1)")
    else:
        print("⚠️  Higher gradients detected - may indicate steep regions")

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)

    # Optional: Save all results for further analysis
    print(f"\nSaving analysis results...")
    import pickle

    results_summary = {
        'maxima': maxima,
        'basin_analysis': basin_analysis,
        'persistent_homology': ph_result,
        'mapper_graph': mapper_result,
        'gradient_analysis': gradient_analysis
    }

    with open('hyperspace_analysis_results.pkl', 'wb') as f:
        pickle.dump(results_summary, f)

    print("Results saved to 'hyperspace_analysis_results.pkl'")
    print("Mapper visualization saved to 'mapper_4d_analysis.png' (if successful)")

