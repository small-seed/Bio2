import sys

class SequenceDotPlot:
    def __init__(self, primary_seq, secondary_seq):
        self.primary = primary_seq
        self.secondary = secondary_seq

    def initialize_matrix(self, cols, rows):
        return [[0 for _ in range(cols)] for _ in range(rows)]

    def generate_dot_matrix(self):
        plot_matrix = self.initialize_matrix(len(self.secondary), len(self.primary))
        for row_idx in range(len(self.primary)):
            for col_idx in range(len(self.secondary)):
                if self.primary[row_idx] == self.secondary[col_idx]:
                    plot_matrix[row_idx][col_idx] = 1
        return plot_matrix

    def apply_window_filter(self, window_size, threshold):
        filtered_matrix = self.initialize_matrix(len(self.secondary), len(self.primary))
        half_window = window_size // 2
        for row in range(half_window, len(self.primary) - half_window):
            for col in range(half_window, len(self.secondary) - half_window):
                match_count = 0
                col_offset = col - half_window
                for row_offset in range(-half_window, half_window + 1):
                    current_row = row + row_offset
                    current_col = col_offset + row_offset
                    if self.primary[current_row] == self.secondary[current_col]:
                        match_count += 1
                if match_count >= threshold:
                    filtered_matrix[row][col] = 1
        return filtered_matrix

    def render_plot(self, matrix):
        sys.stdout.write(" " + self.secondary + "\n")
        for row_idx in range(len(matrix)):
            sys.stdout.write(self.primary[row_idx])
            for col_val in matrix[row_idx]:
                if col_val > 0:
                    sys.stdout.write("*")
                else:
                    sys.stdout.write(" ")
            sys.stdout.write("\n")

def example_usage():
    plotter = SequenceDotPlot("CGATATACCTAG", "TATGATGGATT")
    filtered_result = plotter.apply_window_filter(5, 4)
    plotter.render_plot(filtered_result)

if __name__ == "__main__":
    example_usage()
