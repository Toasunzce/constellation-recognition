from sklearn.neighbors import KDTree 
from itertools import combinations
from collections import Counter
import pandas as pd
import numpy as np
import cv2
from constellation_lines import lines


metrics = pd.read_csv('data/triangles.csv')
tree = KDTree(metrics[['angle1', 'angle2', 'angle3']])


def distance(a: pd.DataFrame, b: pd.DataFrame) -> float:
    return np.sqrt(
        (a['x'].values[0] - b['x'].values[0]) ** 2
         + (a['y'].values[0] - b['y'].values[0]) ** 2
    )


def calculate_2d_triangle_angles(side_ab: float, side_bc: float, side_ac: float) -> np.ndarray:
    cos_A = (side_bc**2 + side_ac**2 - side_ab**2) / (2 * side_bc * side_ac)
    cos_B = (side_ab**2 + side_ac**2 - side_bc**2) / (2 * side_ab * side_ac)
    cos_C = (side_ab**2 + side_bc**2 - side_ac**2) / (2 * side_ab * side_bc)
    
    angles = np.degrees([
        np.arccos(np.clip(cos_A, -1.0, 1.0)),
        np.arccos(np.clip(cos_B, -1.0, 1.0)),
        np.arccos(np.clip(cos_C, -1.0, 1.0))
    ])
    return np.sort(angles)


def find_stars(path: str) -> list[tuple[float, float]]:
    img = cv2.imread(path)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    blurred = cv2.medianBlur(gray, ksize=3)
    _, thresh = cv2.threshold(blurred, 10, 255, cv2.THRESH_BINARY)

    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    star_coords = []

    for contour in contours:
        area = cv2.contourArea(contour)
        if area > 0:
            moments = cv2.moments(contour)
            if moments['m00'] != 0:
                x = int(moments['m10']/moments['m00'])
                y = int(moments['m01']/moments['m00'])
                star_coords.append((x, y))
    return star_coords


def find_lines(matches, lines):
    pairs_to_draw = []
    name_votes = []
    for i, j in combinations(range(len(matches)), 2):
        dot_1 = matches[i]
        dot_2 = matches[j]
        if dot_1 == dot_2:
            continue
        for constellation, pairs in lines.items():
            for pair in pairs:
                if dot_1 in pair and dot_2 in pair:
                    name_votes.append(constellation)
                    pairs_to_draw.append([dots[i], dots[j]])
    name = Counter(name_votes).most_common(1)[0][0]
    return name, pairs_to_draw


def draw_lines(lines_to_draw, name: str) -> None:
    img = cv2.imread(path)
    for pair in lines_to_draw:
        dot_1, dot_2 = pair[0], pair[1]
        cv2.line(img, dot_1, dot_2, (255, 255, 0))
    cv2.imwrite(f'examples/results/{name.lower()}.jpg', img)


def choose_best_hip(idx: int, candidates: list[float],
                    distances: dict[float, list[float]],
                      used_hips: set[float]) -> float | None:
    
    if not candidates:
        return None

    counter = Counter(candidates)
    most_common = counter.most_common()
    top_votes = most_common[0][1]
    top_candidates = [hip for hip, votes in most_common if votes == top_votes and hip not in used_hips]

    if not top_candidates:
        return None

    if len(top_candidates) == 1:
        chosen = top_candidates[0]
        return chosen
    else:
        chosen = min(top_candidates, key=lambda hip: np.mean(distances[hip]))
        return chosen


def match_stars_to_catalogue(df: pd.DataFrame, metrics: pd.DataFrame, tree: KDTree) -> dict[int, float]:
    used_hips: set[float] = set()
    matches: dict[int, float] = {}

    for i in range(len(df)):
        dot_1 = df.iloc[[i]]
        hip_distances: dict[float, list[float]] = {}
        hip_candidates: list[float] = []

        for j, k in combinations(range(len(df)), 2):
            if i in [j, k]:
                continue

            dot_2 = df.iloc[[j]]
            dot_3 = df.iloc[[k]]

            AB = distance(dot_1, dot_2)
            BC = distance(dot_2, dot_3)
            AC = distance(dot_1, dot_3)

            angles = calculate_2d_triangle_angles(AB, BC, AC)
            dist, ind = tree.query([angles], k=2)

            for idx, d in zip(ind[0], dist[0]):
                triangle = metrics.iloc[idx]
                hips = [triangle['hip1'], triangle['hip2'], triangle['hip3']]
                for hip in hips:
                    if hip in used_hips:
                        continue
                    hip_candidates.append(hip)
                    hip_distances.setdefault(hip, []).append(d)

        best_hip = choose_best_hip(i, hip_candidates, hip_distances, used_hips)
        if best_hip is not None:
            matches[i] = float(best_hip)
            used_hips.add(best_hip)

    return matches




path = 'examples/inputs/image1.png'

dots = find_stars(path)
df = pd.DataFrame(dots, columns=["x", "y"])
matches = match_stars_to_catalogue(df, metrics, tree)

name, output_list = find_lines(matches, lines)
draw_lines(output_list, name)



# TODO: Обрабатывать случаи с несколькими созвездиями на одном фото.
# Идея: сопоставлять треугольники между всеми звездами, а не только внутри одного созвездия.
