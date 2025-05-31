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
# import kmapper
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


class ReactionHyperspaceAnalyzer:
    """
    Geometric analysis of reaction condition spaces

    Designed for experimental datasets where reaction conditions (temperature,
    concentrations, etc.) map to reaction yields. Handles both 3D and 4D
    hyperspace analysis with flexible column naming.
    """

    def __init__(self, csv_path: str, coordinate_columns: List[str], yield_column: str, coordinate_column_labels=None):
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
        if coordinate_column_labels is not None:
            self.coordinate_column_labels = coordinate_column_labels
        else:
            self.coordinate_column_labels = self.coordinate_columns
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
        Dict containing:
            'gradient_magnitude': N-D array of gradient magnitudes
            'gradient_components': List of N-D arrays (one per dimension)
            'mean_gradient': Average gradient magnitude
            'max_gradient': Maximum gradient magnitude
        """
        if self.yield_array is None:
            raise ValueError("Must call interpolate_to_regular_grid() first")

        print("Computing yield gradients...")

        # Compute gradients using numpy.gradient (central differences)
        gradients = np.gradient(self.yield_array)

        # Convert to real coordinate units (not grid units)
        grid_spacings = []
        for i in range(self.dimensionality):
            coord_range = self.coordinate_ranges[i][1] - self.coordinate_ranges[i][0]
            grid_spacing = coord_range / (self.grid_resolution - 1)
            grid_spacings.append(grid_spacing)
            gradients[i] = gradients[i] / grid_spacing

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

            # Compute robustness as average yield in basin
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

    def analyze_transition_paths(self, maxima_result, n_path_points=100, method='comprehensive'):
        """
        Find and analyze transition paths between identified maxima.

        Addresses reviewer requests for:
        - Systematic algorithms for hyperspace structure analysis
        - Geometric insights into high-dimensional spaces
        - Characterization of disconnected maxima regions

        Parameters:
        -----------
        maxima_result : Dict
            Output from find_local_maxima()
        n_path_points : int, default=100
            Number of points to sample along transition paths
        method : str, default='comprehensive'
            Analysis method ('direct', 'saddle_search', 'comprehensive')

        Returns:
        --------
        Dict containing transition path analysis and visualizations
        """
        if self.yield_array is None:
            raise ValueError("Must call interpolate_to_regular_grid() first")

        if maxima_result['count'] < 2:
            print("Need at least 2 maxima for transition path analysis")
            return None

        print(f"Analyzing transition paths between {maxima_result['count']} maxima...")

        # Extract peak information
        peak1_coords = maxima_result['coordinates'][0]
        peak2_coords = maxima_result['coordinates'][1]
        peak1_yield = maxima_result['yields'][0]
        peak2_yield = maxima_result['yields'][1]
        peak1_idx = maxima_result['indices'][0]
        peak2_idx = maxima_result['indices'][1]

        print(f"Peak 1: yield={peak1_yield:.4f} at {dict(zip(self.coordinate_columns, peak1_coords))}")
        print(f"Peak 2: yield={peak2_yield:.4f} at {dict(zip(self.coordinate_columns, peak2_coords))}")

        # 1. Direct linear path analysis
        print("Analyzing direct linear path...")
        path_analysis = self._analyze_direct_path(peak1_coords, peak2_coords, peak1_idx, peak2_idx, n_path_points)

        # 2. Find saddle points and barriers
        print("Searching for saddle points and barriers...")
        saddle_analysis = self._find_transition_barriers(peak1_idx, peak2_idx, path_analysis)

        # 3. Compute geometric properties
        print("Computing path geometry...")
        geometry_analysis = self._compute_path_geometry(path_analysis, peak1_coords, peak2_coords)

        # 4. Create comprehensive visualizations
        print("Creating transition path visualizations...")
        self._visualize_transition_paths(path_analysis, saddle_analysis, geometry_analysis,
                                         peak1_coords, peak2_coords, peak1_yield, peak2_yield)

        # 5. Generate summary report
        summary = self._generate_transition_summary(path_analysis, saddle_analysis, geometry_analysis)

        print("\nTransition Path Analysis Summary:")
        print(f"  Direct path length: {geometry_analysis['path_length_real']:.6f} units")
        print(f"  Barrier height: {saddle_analysis['barrier_height']:.4f} yield units")
        print(f"  Barrier location: {saddle_analysis['barrier_fraction']:.2f} along path")
        print(f"  Path curvature: {geometry_analysis['mean_curvature']:.4f}")
        print(f"  Transition efficiency: {summary['transition_efficiency']:.3f}")

        return {
            'path_analysis': path_analysis,
            'saddle_analysis': saddle_analysis,
            'geometry_analysis': geometry_analysis,
            'summary': summary,
            'peak1_coords': peak1_coords,
            'peak2_coords': peak2_coords,
            'peak1_yield': peak1_yield,
            'peak2_yield': peak2_yield
        }

    def _analyze_direct_path(self, peak1_coords, peak2_coords, peak1_idx, peak2_idx, n_points):
        """Analyze direct linear path between two peaks."""

        # Create linear path in parameter space
        path_coords = np.zeros((n_points, self.dimensionality))
        path_yields = np.zeros(n_points)
        path_indices = np.zeros((n_points, self.dimensionality), dtype=int)

        for i in range(n_points):
            t = i / (n_points - 1)  # Parameter from 0 to 1

            # Linear interpolation in parameter space
            current_coords = (1 - t) * peak1_coords + t * peak2_coords
            path_coords[i] = current_coords

            # Convert to grid indices for yield lookup
            grid_indices = np.zeros(self.dimensionality, dtype=int)
            for dim in range(self.dimensionality):
                # Map coordinate to grid index
                coord_min, coord_max = self.coordinate_ranges[dim]
                normalized = (current_coords[dim] - coord_min) / (coord_max - coord_min)
                grid_idx = int(np.clip(normalized * (self.grid_resolution - 1), 0, self.grid_resolution - 1))
                grid_indices[dim] = grid_idx

            path_indices[i] = grid_indices

            # Look up yield at this point
            path_yields[i] = self.yield_array[tuple(grid_indices)]

        return {
            'path_coords': path_coords,
            'path_yields': path_yields,
            'path_indices': path_indices,
            'n_points': n_points,
            'path_parameters': np.linspace(0, 1, n_points)
        }

    def _find_transition_barriers(self, peak1_idx, peak2_idx, path_analysis):
        """Find barriers and saddle-like regions along the transition path."""

        path_yields = path_analysis['path_yields']
        path_coords = path_analysis['path_coords']
        path_params = path_analysis['path_parameters']

        # Find minimum yield point along path (main barrier)
        barrier_idx = np.argmin(path_yields)
        barrier_yield = path_yields[barrier_idx]
        barrier_coords = path_coords[barrier_idx]
        barrier_param = path_params[barrier_idx]

        # Compute barrier height relative to both peaks
        peak1_yield = path_yields[0]
        peak2_yield = path_yields[-1]
        barrier_height = min(peak1_yield, peak2_yield) - barrier_yield

        # Find gradient at barrier point
        barrier_gradient = self._compute_local_gradient(path_analysis['path_indices'][barrier_idx])
        gradient_magnitude = np.linalg.norm(barrier_gradient)

        # Analyze barrier width (region within 90% of barrier depth)
        barrier_threshold = barrier_yield + 0.1 * barrier_height
        barrier_region = path_yields <= barrier_threshold
        barrier_width_indices = np.where(barrier_region)[0]
        barrier_width = len(barrier_width_indices) / len(path_yields) if len(barrier_width_indices) > 0 else 0

        return {
            'barrier_idx': barrier_idx,
            'barrier_yield': barrier_yield,
            'barrier_coords': barrier_coords,
            'barrier_fraction': barrier_param,
            'barrier_height': barrier_height,
            'barrier_gradient': barrier_gradient,
            'gradient_magnitude': gradient_magnitude,
            'barrier_width': barrier_width,
            'is_saddle_like': gradient_magnitude < 0.01  # Low gradient suggests saddle-like point
        }

    def _compute_local_gradient(self, grid_indices):
        """Compute gradient at a specific grid point using finite differences."""

        gradient = np.zeros(self.dimensionality)

        for dim in range(self.dimensionality):
            # Forward and backward indices
            idx_forward = grid_indices.copy()
            idx_backward = grid_indices.copy()

            # Finite difference with boundary checking
            if grid_indices[dim] < self.grid_resolution - 1:
                idx_forward[dim] += 1
                forward_yield = self.yield_array[tuple(idx_forward)]
            else:
                forward_yield = self.yield_array[tuple(grid_indices)]

            if grid_indices[dim] > 0:
                idx_backward[dim] -= 1
                backward_yield = self.yield_array[tuple(idx_backward)]
            else:
                backward_yield = self.yield_array[tuple(grid_indices)]

            # Convert to real coordinate spacing
            coord_range = self.coordinate_ranges[dim][1] - self.coordinate_ranges[dim][0]
            spacing = coord_range / (self.grid_resolution - 1)

            # Central difference approximation
            gradient[dim] = (forward_yield - backward_yield) / (2 * spacing)

        return gradient

    def _compute_path_geometry(self, path_analysis, peak1_coords, peak2_coords):
        """Compute geometric properties of the transition path."""

        path_coords = path_analysis['path_coords']
        path_yields = path_analysis['path_yields']

        # Compute path length in parameter space
        path_segments = np.diff(path_coords, axis=0)
        segment_lengths = np.linalg.norm(path_segments, axis=1)
        path_length_real = np.sum(segment_lengths)

        # Compute path curvature (change in direction)
        if len(path_coords) >= 3:
            # Second derivatives for curvature
            second_derivatives = np.diff(path_coords, n=2, axis=0)
            curvatures = np.linalg.norm(second_derivatives, axis=1)
            mean_curvature = np.mean(curvatures)
            max_curvature = np.max(curvatures)
        else:
            mean_curvature = max_curvature = 0

        # Compute yield gradient along path
        yield_gradient = np.gradient(path_yields)
        max_yield_gradient = np.max(np.abs(yield_gradient))

        # Direct distance between peaks
        direct_distance = np.linalg.norm(peak2_coords - peak1_coords)
        path_efficiency = direct_distance / path_length_real if path_length_real > 0 else 0

        return {
            'path_length_real': path_length_real,
            'direct_distance': direct_distance,
            'path_efficiency': path_efficiency,
            'mean_curvature': mean_curvature,
            'max_curvature': max_curvature,
            'max_yield_gradient': max_yield_gradient,
            'segment_lengths': segment_lengths
        }

    def _generate_transition_summary(self, path_analysis, saddle_analysis, geometry_analysis):
        """Generate summary metrics for transition analysis."""

        path_yields = path_analysis['path_yields']

        # Transition efficiency metrics
        yield_range = np.max(path_yields) - np.min(path_yields)
        barrier_prominence = saddle_analysis['barrier_height'] / yield_range if yield_range > 0 else 0
        transition_efficiency = 1.0 - barrier_prominence  # Higher when barrier is lower

        # Geometric complexity
        geometric_complexity = geometry_analysis['mean_curvature'] * geometry_analysis['path_length_real']

        # Classification
        if saddle_analysis['is_saddle_like']:
            transition_type = "Saddle-mediated transition"
        elif saddle_analysis['barrier_height'] < 0.01:
            transition_type = "Nearly barrierless transition"
        else:
            transition_type = "Barrier-mediated transition"

        return {
            'transition_efficiency': transition_efficiency,
            'barrier_prominence': barrier_prominence,
            'geometric_complexity': geometric_complexity,
            'transition_type': transition_type,
            'yield_range_along_path': yield_range
        }

    def _visualize_transition_paths(self, path_analysis, saddle_analysis, geometry_analysis,
                                    peak1_coords, peak2_coords, peak1_yield, peak2_yield):
        """Create comprehensive visualizations of transition paths."""

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        # Create figure with multiple subplots
        fig = plt.figure( figsize=(15, 15))

        # 1. 3D path trajectory (use first 3 dimensions for visualization)
        ax1 = fig.add_subplot(2, 2, 1, projection='3d')
        self._plot_3d_trajectory(ax1, path_analysis, peak1_coords, peak2_coords,
                                 peak1_yield, peak2_yield, saddle_analysis)

        # 2. Yield profile along path
        ax2 = fig.add_subplot(2, 2, 2)
        self._plot_yield_profile(ax2, path_analysis, saddle_analysis)

        # 3. Parameter evolution along path
        ax3 = fig.add_subplot(2, 2, 3)
        self._plot_parameter_evolution(ax3, path_analysis)

        # 4. 2D cross-sections through barrier region
        ax4 = fig.add_subplot(2, 2, 4)
        self._plot_barrier_cross_section(ax4, path_analysis, saddle_analysis)

        # # 5. Path geometry analysis
        # ax5 = fig.add_subplot(2, 3, 5)
        # self._plot_path_geometry(ax5, path_analysis, geometry_analysis)
        #
        # # 6. Summary metrics visualization
        # ax6 = fig.add_subplot(2, 3, 6)
        # self._plot_summary_metrics(ax6, path_analysis, saddle_analysis, geometry_analysis)

        # plt.suptitle('Transition Path Analysis Between Reaction Maxima', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig('transition_path_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()

    def _plot_3d_trajectory(self, ax, path_analysis, peak1_coords, peak2_coords,
                            peak1_yield, peak2_yield, saddle_analysis):
        """Plot 3D trajectory of transition path using first 3 coordinate dimensions."""

        path_coords = path_analysis['path_coords']
        path_yields = path_analysis['path_yields']

        # Use first 3 dimensions for 3D visualization
        x = path_coords[:, 0]
        y = path_coords[:, 1]
        z = path_coords[:, 2]

        # Color by yield
        scatter = ax.scatter(x, y, z, c=path_yields, cmap='viridis', s=50, alpha=0.7)

        # Plot path as line
        ax.plot(x, y, z, 'k-', alpha=0.3, linewidth=2)

        # Highlight peaks
        ax.scatter([peak1_coords[0]], [peak1_coords[1]], [peak1_coords[2]],
                   c='green', s=200, marker='*', label=f'Global maximum (yield={peak1_yield:.3f})', edgecolors='black')
        ax.scatter([peak2_coords[0]], [peak2_coords[1]], [peak2_coords[2]],
                   c='blue', s=200, marker='*', label=f'Local maximum (yield={peak2_yield:.3f})', edgecolors='black')

        # Highlight barrier point
        barrier_idx = saddle_analysis['barrier_idx']
        barrier_coords = saddle_analysis['barrier_coords']
        ax.scatter([barrier_coords[0]], [barrier_coords[1]], [barrier_coords[2]],
                   c='red', s=150, marker='X', label=f'Barrier (yield={saddle_analysis["barrier_yield"]:.3f})',
                   edgecolors='black', zorder=10)

        # Labels and formatting
        ax.set_xlabel(self.coordinate_column_labels[0])
        ax.set_ylabel(self.coordinate_column_labels[1])
        ax.set_zlabel(self.coordinate_column_labels[2])
        ax.set_title('3D projection of 4D transition trajectory')
        ax.legend()

        # Add colorbar
        plt.colorbar(scatter, ax=ax, shrink=0.8, label='Yield')

    def _plot_yield_profile(self, ax, path_analysis, saddle_analysis):
        """Plot yield profile along the transition path."""

        path_params = path_analysis['path_parameters']
        path_yields = path_analysis['path_yields']

        # Plot yield profile
        ax.plot(path_params, path_yields, 'b-', linewidth=3, label='Yield along path')

        # Highlight barrier point
        barrier_param = saddle_analysis['barrier_fraction']
        barrier_yield = saddle_analysis['barrier_yield']
        ax.plot(barrier_param, barrier_yield, 'rx', label='Barrier point', markersize=15, markeredgewidth=3)

        # Highlight peaks
        ax.plot(0, path_yields[0], '*g', markersize=15, label='Global maximum')
        ax.plot(1, path_yields[-1], '*b', markersize=15, label='Local maximum')

        # Add barrier height annotation
        barrier_height = saddle_analysis['barrier_height']
        ax.annotate(f'Barrier\n(point of lowest yield, {barrier_height:.4f})',
                    xy=(barrier_param, barrier_yield), xytext=(0.5, barrier_yield + 0.05),
                    arrowprops=dict(arrowstyle='->', color='red'),
                    fontsize=10, ha='center')

        ax.set_xlabel('Path internal coordinate (0=Peak1, 1=Peak2)')
        ax.set_ylabel('Yield')
        ax.set_title('Yield profile along transition path')
        ax.legend()
        ax.grid(True, alpha=0.3)

    def _plot_parameter_evolution(self, ax, path_analysis):
        """Plot how each parameter changes along the transition path."""

        path_params = path_analysis['path_parameters']
        path_coords = path_analysis['path_coords']

        # Plot each coordinate dimension
        colors = ['red', 'blue', 'green', 'purple', 'orange']
        for dim in range(self.dimensionality):
            ax.plot(path_params, path_coords[:, dim],
                    color=colors[dim % len(colors)], linewidth=2,
                    label=self.coordinate_column_labels[dim])

        ax.set_xlabel('Path internal coordinate (0=Global max., 1=local max.)')
        ax.set_ylabel('Parameter value')
        ax.set_title('Parameter evolution along transition path')
        ax.legend()
        ax.grid(True, alpha=0.3)

    def _plot_barrier_cross_section(self, ax, path_analysis, saddle_analysis):
        """Plot 2D cross-section around the barrier region."""

        barrier_coords = saddle_analysis['barrier_coords']

        # Create 2D slice through the barrier point using first 2 dimensions
        n_points = 50
        dim1, dim2 = 0, 1  # Use first two dimensions

        # Create grid around barrier point
        coord_ranges = []
        for dim in [dim1, dim2]:
            coord_min, coord_max = self.coordinate_ranges[dim]
            center = barrier_coords[dim]
            half_width = (coord_max - coord_min) * 0.3  # 30% of total range
            grid_min = max(coord_min, center - half_width)
            grid_max = min(coord_max, center + half_width)
            coord_ranges.append(np.linspace(grid_min, grid_max, n_points))

        X, Y = np.meshgrid(coord_ranges[0], coord_ranges[1])
        Z = np.zeros_like(X)

        # Sample yields on this 2D grid (keeping other dimensions at barrier values)
        for i in range(n_points):
            for j in range(n_points):
                sample_coords = barrier_coords.copy()
                sample_coords[dim1] = X[i, j]
                sample_coords[dim2] = Y[i, j]

                # Convert to grid indices
                grid_indices = []
                for dim in range(self.dimensionality):
                    coord_min, coord_max = self.coordinate_ranges[dim]
                    normalized = (sample_coords[dim] - coord_min) / (coord_max - coord_min)
                    grid_idx = int(np.clip(normalized * (self.grid_resolution - 1), 0, self.grid_resolution - 1))
                    grid_indices.append(grid_idx)

                Z[i, j] = self.yield_array[tuple(grid_indices)]

        # Create contour plot
        contour = ax.contourf(X, Y, Z, levels=20, cmap='viridis', alpha=0.8, method='cubic')
        ax.contour(X, Y, Z, levels=10, colors='black', alpha=0.3, linewidths=0.5, method='cubic')

        # Mark barrier point
        ax.plot(barrier_coords[dim1], barrier_coords[dim2], 'rx', markersize=15, markeredgewidth=3)

        ax.set_xlabel(self.coordinate_column_labels[dim1])
        ax.set_ylabel(self.coordinate_column_labels[dim2])
        ax.set_title(
            f'Cross-section at barrier point\n({self.coordinate_column_labels[dim1]} vs {self.coordinate_column_labels[dim2]})')

        plt.colorbar(contour, ax=ax, label='Yield')

    def _plot_path_geometry(self, ax, path_analysis, geometry_analysis):
        """Plot geometric properties of the path."""

        path_params = path_analysis['path_parameters']
        segment_lengths = geometry_analysis['segment_lengths']

        # Plot cumulative distance along path
        cumulative_distance = np.cumsum(np.concatenate([[0], segment_lengths]))
        ax.plot(path_params, cumulative_distance, 'g-', linewidth=2, label='Cumulative distance')

        # Plot local curvature if available
        if len(path_analysis['path_coords']) >= 3:
            path_coords = path_analysis['path_coords']
            second_derivatives = np.diff(path_coords, n=2, axis=0)
            curvatures = np.linalg.norm(second_derivatives, axis=1)
            curvature_params = path_params[1:-1]  # Curvature has 2 fewer points

            ax2 = ax.twinx()
            ax2.plot(curvature_params, curvatures, 'r--', linewidth=2, label='Path curvature')
            ax2.set_ylabel('Curvature', color='red')
            ax2.tick_params(axis='y', labelcolor='red')

        ax.set_xlabel('Path Parameter')
        ax.set_ylabel('Cumulative Distance', color='green')
        ax.set_title('Path Geometry Analysis')
        ax.tick_params(axis='y', labelcolor='green')
        ax.grid(True, alpha=0.3)

    def _plot_summary_metrics(self, ax, path_analysis, saddle_analysis, geometry_analysis):
        """Plot summary metrics as a text-based visualization."""

        ax.axis('off')  # Turn off axes for text display

        # Compile summary statistics
        summary_text = f"""
        TRANSITION PATH ANALYSIS SUMMARY
        
        Path Properties:
        • Total path length: {geometry_analysis['path_length_real']:.6f} units
        • Direct distance: {geometry_analysis['direct_distance']:.6f} units  
        • Path efficiency: {geometry_analysis['path_efficiency']:.3f}
        • Mean curvature: {geometry_analysis['mean_curvature']:.4f}
        
        Barrier Analysis:
        • Barrier height: {saddle_analysis['barrier_height']:.4f} yield units
        • Barrier position: {saddle_analysis['barrier_fraction']:.2f} along path
        • Barrier yield: {saddle_analysis['barrier_yield']:.4f}
        • Gradient magnitude: {saddle_analysis['gradient_magnitude']:.4f}
        • Saddle-like: {'Yes' if saddle_analysis['is_saddle_like'] else 'No'}
        
        Peak Information:
        • Peak 1 coordinates: {dict(zip(self.coordinate_columns, path_analysis['path_coords'][0]))}
        • Peak 2 coordinates: {dict(zip(self.coordinate_columns, path_analysis['path_coords'][-1]))}
        
        Geometric Classification:
        • Transition type: {'Saddle-mediated' if saddle_analysis['is_saddle_like'] else 'Barrier-mediated'}
        • Yield range: {np.max(path_analysis['path_yields']) - np.min(path_analysis['path_yields']):.4f}
        """

        ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))


# Usage example for testing
if __name__ == "__main__":
    pass
    # For 4D analysis
    analyzer = ReactionHyperspaceAnalyzer(
        csv_path='postprocessed_yields_decimated.csv',
        coordinate_columns=['ic001','am001','ald001','ptsa'],
        yield_column='yield',
        coordinate_column_labels=['Isocyanide', 'Amine', 'Aldehyde', 'pTSA'],
    )

    # Interpolate to regular grid
    grid_resolution = 40
    analyzer.interpolate_to_regular_grid(grid_resolution=grid_resolution)

    analyzer.save_regular_grid_to_file('4d_grid.npz')

    analyzer.load_regular_grid_from_file('4d_grid.npz')

    # Validate quality
    analyzer.validate_grid_quality()

    # Find local maxima
    maxima = analyzer.find_local_maxima(min_distance=int(round(grid_resolution/2)), threshold_abs=0.11)

    # # Compute gradients
    # gradients = analyzer.compute_yield_gradients()

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

    # ================== TRANSITION ANALYSIS ==================
    analyzer = ReactionHyperspaceAnalyzer(
        csv_path='postprocessed_yields_decimated.csv',
        coordinate_columns=['ic001','am001','ald001','ptsa'],
        yield_column='yield',
        coordinate_column_labels=['Isocyanide', 'Amine', 'Aldehyde', 'pTSA'],
    )

    # Interpolate to regular grid
    analyzer.interpolate_to_regular_grid(grid_resolution=19)

    analyzer.save_regular_grid_to_file('4d_grid.npz')

    analyzer.load_regular_grid_from_file('4d_grid.npz')

    # Validate quality
    analyzer.validate_grid_quality()

    # Find local maxima
    maxima = analyzer.find_local_maxima(min_distance=10, threshold_abs=0.11)


    # After your existing basin analysis:
    transition_analysis = analyzer.analyze_transition_paths(maxima, n_path_points=50)
