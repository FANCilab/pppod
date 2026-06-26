from pathlib import Path
import os
import re
from datetime import datetime
import numpy as np


from pathlib import Path
import re


from pathlib import Path
from datetime import datetime
import re
import numpy as np
from dateutil import parser as date_parser

def s2p_run_paths(root_folder):
    """
    Find Suite2p output folders ready for ALF conversion.

    A Suite2p output folder is defined as the parent folder of one or more
    plane folders named plane0, plane1, plane2, etc.

    A Suite2p output folder is considered ready only if:
      1. A .log file exists alongside the planeX folders.
      2. db.npy exists alongside the planeX folders.
      3. Every planeX folder contains:
            F.npy
            Fneu.npy
            spks.npy
            stat.npy
            iscell.npy

    Parameters
    ----------
    root_folder : str or pathlib.Path
        Root folder to recursively scan.

    Returns
    -------
    list[pathlib.Path]
        List of Suite2p parent folders ready for ALF conversion.
    """

    root_folder = Path(root_folder)

    required_files = {
        "F.npy",
        "Fneu.npy",
        "stat.npy",
        "spks.npy",
        "iscell.npy",
    }

    plane_pattern = re.compile(r"^plane\d+$", re.IGNORECASE)

    if not root_folder.exists():
        raise FileNotFoundError(f"Root folder does not exist: {root_folder}")

    # Find all planeX folders under root_folder, including if root_folder itself
    # is the Suite2p output folder.
    plane_folders = [
        p for p in root_folder.rglob("*")
        if p.is_dir() and plane_pattern.match(p.name)
    ]

    # If root_folder itself directly contains plane folders, rglob("*") catches them.
    # Group plane folders by their parent Suite2p run folder.
    candidate_parents = {}

    for plane_folder in plane_folders:
        parent = plane_folder.parent

        # Ignore combined/planeX if such a strange thing exists.
        if parent.name.lower() == "combined":
            continue

        candidate_parents.setdefault(parent, []).append(plane_folder)

    ready_suite2p_folders = []

    for parent, planes in candidate_parents.items():
        # Must have db.npy alongside plane0, plane1, etc.
        if not (parent / "db.npy").is_file():
            continue

        # Must have run.log alongside plane0, plane1, etc.
        if not (parent / "run.log").is_file():
            continue


        all_planes_ready = True

        for plane in planes:
            try:
                files_in_plane = {
                    p.name for p in plane.iterdir()
                    if p.is_file()
                }
            except OSError:
                all_planes_ready = False
                break

            if not required_files.issubset(files_in_plane):
                all_planes_ready = False
                break

        if all_planes_ready:
            ready_suite2p_folders.append(parent)

    return sorted(ready_suite2p_folders)

def s2p_filter_run_paths(s2p_run_paths):
    """
    DEPRICATED FUNCTION DO NOT USE AT ANY COST
    """

    from pathlib import Path
    from datetime import datetime
    import re
    import numpy as np
    from dateutil import parser as date_parser

    valid_paths = []
    session_info = []
    data_path_lists = []

    for run_path in s2p_run_paths:
        run_path = Path(run_path)

        # ------------------------------------------------------------
        # 1. Check that the Suite2p run path itself contains exactly one
        #    date-like folder component.
        # ------------------------------------------------------------
        run_parts = [
            p for p in re.split(r"[\\/]+", str(run_path))
            if p not in ("", ".", "..")
        ]

        run_dates = []

        for part in run_parts:
            token = part.strip()

            has_4digit_year = re.search(r"(19|20)\d{2}", token) is not None
            has_compact_yyyymmdd = re.fullmatch(r"(19|20)\d{6}", token) is not None
            has_month_name = re.search(
                r"(jan|feb|mar|apr|may|jun|jul|aug|sep|sept|oct|nov|dec)",
                token,
                re.IGNORECASE,
            ) is not None

            if not (has_4digit_year or has_compact_yyyymmdd or has_month_name):
                continue

            try:
                parsed_date = date_parser.parse(
                    token,
                    fuzzy=False,
                    default=datetime(1900, 1, 1),
                )
            except Exception:
                continue

            if parsed_date.year == 1900:
                continue

            run_dates.append(parsed_date.strftime("%Y-%m-%d"))

        run_dates = sorted(set(run_dates))

        if len(run_dates) != 1:
            print(
                f"[EXCLUDED] {run_path} | "
                f"Suite2p run path must contain exactly one date-like folder; "
                f"found {len(run_dates)}: {run_dates}"
            )
            continue

        run_date = run_dates[0]

        # ------------------------------------------------------------
        # 2. Load top-level db.npy. This must be alongside planeX folders,
        #    not inside planeX.
        # ------------------------------------------------------------
        db_path = run_path / "db.npy"

        if not db_path.is_file():
            print(
                f"[EXCLUDED] {run_path} | "
                f"Missing top-level db.npy alongside planeX folders: {db_path}"
            )
            continue

        try:
            db = np.load(db_path, allow_pickle=True).item()
        except Exception as exc:
            print(
                f"[EXCLUDED] {run_path} | "
                f"Could not load db.npy: {db_path}; error: {exc}"
            )
            continue

        if not isinstance(db, dict):
            print(
                f"[EXCLUDED] {run_path} | "
                f"db.npy did not load as a dictionary; got type {type(db)}"
            )
            continue

        data_paths = db.get("data_path", None)

        if data_paths is None:
            print(
                f"[EXCLUDED] {run_path} | "
                f'db.npy does not contain required key "data_path"'
            )
            continue

        if isinstance(data_paths, (str, Path)):
            data_paths = [data_paths]
        else:
            try:
                data_paths = list(data_paths)
            except TypeError:
                data_paths = [data_paths]

        data_paths = [str(p) for p in data_paths]

        if len(data_paths) == 0:
            print(
                f"[EXCLUDED] {run_path} | "
                f'db["data_path"] exists but is empty'
            )
            continue

        extracted_subject_date_session = []
        exclusion_reason = None

        # ------------------------------------------------------------
        # 3. Every db["data_path"] entry must resolve to:
        #       ... / Subject / Date / Session
        # ------------------------------------------------------------
        for raw_path_string in data_paths:
            raw_path_string = str(raw_path_string).strip()

            if not raw_path_string:
                extracted_subject_date_session = []
                exclusion_reason = 'db["data_path"] contains an empty path entry'
                break

            raw_parts = [
                p for p in re.split(r"[\\/]+", raw_path_string)
                if p not in ("", ".", "..")
            ]

            date_hits = []

            for i, part in enumerate(raw_parts):
                token = part.strip()

                has_4digit_year = re.search(r"(19|20)\d{2}", token) is not None
                has_compact_yyyymmdd = re.fullmatch(r"(19|20)\d{6}", token) is not None
                has_month_name = re.search(
                    r"(jan|feb|mar|apr|may|jun|jul|aug|sep|sept|oct|nov|dec)",
                    token,
                    re.IGNORECASE,
                ) is not None

                if not (has_4digit_year or has_compact_yyyymmdd or has_month_name):
                    continue

                try:
                    parsed_date = date_parser.parse(
                        token,
                        fuzzy=False,
                        default=datetime(1900, 1, 1),
                    )
                except Exception:
                    continue

                if parsed_date.year == 1900:
                    continue

                formatted_date = parsed_date.strftime("%Y-%m-%d")
                date_hits.append((i, formatted_date))

            # Each raw data path usually contains exactly one date folder. Some
            # db.py files store a shortened .../2P/subject/session path; in that
            # case use the already-validated Suite2p run-path date.
            if len(date_hits) == 0:
                two_p_indices = [
                    i for i, part in enumerate(raw_parts)
                    if part.lower() == "2p"
                ]
                if not two_p_indices:
                    extracted_subject_date_session = []
                    exclusion_reason = (
                        f'Could not infer missing date for db["data_path"] entry; '
                        f'path "{raw_path_string}" had no date-like folder and no "2P" folder'
                    )
                    break

                two_p_index = two_p_indices[-1]
                if two_p_index + 2 >= len(raw_parts):
                    extracted_subject_date_session = []
                    exclusion_reason = (
                        f'Could not infer subject/session for date-less db["data_path"] entry; '
                        f'expected .../2P/Subject/Session in "{raw_path_string}"'
                    )
                    break

                subject = raw_parts[two_p_index + 1]
                formatted_date = run_date
                session = raw_parts[two_p_index + 2]
            elif len(date_hits) != 1:
                extracted_subject_date_session = []
                exclusion_reason = (
                    f'Each db["data_path"] entry must contain exactly one date-like folder; '
                    f'path "{raw_path_string}" had {len(date_hits)} date hits: {date_hits}'
                )
                break
            else:
                date_index, formatted_date = date_hits[0]

                # Raw data path must be:
                #   ... / Subject / Date / Session
                if date_index == 0 or date_index + 1 >= len(raw_parts):
                    extracted_subject_date_session = []
                    exclusion_reason = (
                        f'Raw data path must be formatted as .../Subject/Date/Session; '
                        f'could not extract subject/session from "{raw_path_string}"'
                    )
                    break

                subject = raw_parts[date_index - 1]
                session = raw_parts[date_index + 1]

            if not subject or not session:
                extracted_subject_date_session = []
                exclusion_reason = (
                    f'Extracted empty subject or session from db["data_path"] entry: '
                    f'"{raw_path_string}"'
                )
                break

            extracted_subject_date_session.append(
                (subject, formatted_date, session)
            )

        if not extracted_subject_date_session:
            if exclusion_reason is None:
                exclusion_reason = (
                    f'Could not extract any valid subject/date/session tuple '
                    f'from db["data_path"]: {data_paths}'
                )

            print(f"[EXCLUDED] {run_path} | {exclusion_reason}")
            continue

        subjects = {x[0] for x in extracted_subject_date_session}
        dates = {x[1] for x in extracted_subject_date_session}
        sessions = {x[2] for x in extracted_subject_date_session}

        # All db["data_path"] entries must agree on subject/date. The session
        # order is preserved because Suite2p concatenates in db["data_path"]
        # order.
        if len(subjects) != 1 or len(dates) != 1:
            print(
                f"[EXCLUDED] {run_path} | "
                f'db["data_path"] entries do not agree on one subject/date; '
                f"subjects={sorted(subjects)}, dates={sorted(dates)}, sessions={sorted(sessions)}"
            )
            continue

        subject = next(iter(subjects))
        date = next(iter(dates))
        session = "_".join(raw_session for _, _, raw_session in extracted_subject_date_session)

        # The date in the Suite2p output path must match the date in db.npy.
        if date != run_date:
            print(
                f"[EXCLUDED] {run_path} | "
                f"Date mismatch between Suite2p run path and db.npy data_path; "
                f"run path date={run_date}, db data_path date={date}"
            )
            continue

        trimmed_data_paths = [
            f"{raw_subject}/{raw_date}/{raw_session}"
            for raw_subject, raw_date, raw_session in extracted_subject_date_session
        ]

        valid_paths.append(run_path)
        session_info.append((subject, date, session))
        data_path_lists.append(trimmed_data_paths)

    return valid_paths, session_info, data_path_lists
