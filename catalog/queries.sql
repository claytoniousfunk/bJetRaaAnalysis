-- Example queries against scan_catalog.db
-- Run with: sqlite3 catalog/scan_catalog.db < catalog/queries.sql
-- or open interactively: sqlite3 catalog/scan_catalog.db


SELECT path, produced_on
  FROM (
      SELECT path, config_key, produced_on,
             ROW_NUMBER() OVER (PARTITION BY config_key ORDER BY produced_on DESC) AS rn
      FROM scan_files
      WHERE dataset_group = 'PbPb' OR dataset_group = 'pp'
  )
  WHERE rn = 1;




-- 1. Overview: file count, total size, date range per dataset group.
-- SELECT
--     dataset_group,
--     COUNT(*)                      AS n_files,
--     ROUND(SUM(size_bytes) / 1e6, 1) AS total_mb,
--     MIN(produced_on)              AS earliest,
--     MAX(produced_on)              AS latest
-- FROM scan_files
-- GROUP BY dataset_group
-- ORDER BY dataset_group;

-- 2. Duplicate configs: same cuts/tags, produced on more than one date.
-- (config_key is the filename with the date token stripped out)
-- SELECT
--     config_key,
--     COUNT(*)                       AS n_versions,
--     GROUP_CONCAT(produced_on, ', ') AS dates
-- FROM scan_files
-- GROUP BY config_key
-- HAVING COUNT(*) > 1;

-- 3. Latest file per config (window function) -- the one you should
--    actually be pointing your plotting macros at.
-- SELECT path, config_key, produced_on
-- FROM (
--     SELECT
--         path, config_key, produced_on,
--         ROW_NUMBER() OVER (PARTITION BY config_key ORDER BY produced_on DESC) AS rn
--     FROM scan_files
-- )
-- WHERE rn = 1;

-- 4. Selection efficiency (events after / events before) per file, as a
--    quick QA ranking -- flags scans with unexpectedly low survival rates.
-- SELECT
--     filename,
--     events_before_selection,
--     events_after_selection,
--     ROUND(events_after_selection / NULLIF(events_before_selection, 0), 4) AS efficiency
-- FROM scan_files
-- WHERE events_before_selection IS NOT NULL
-- ORDER BY efficiency ASC;

-- 5. Tag search: files that have ALL of a given set of cut tags applied
--    (classic "match every row in a set" pattern via GROUP BY/HAVING).
-- SELECT sf.path
-- FROM scan_files sf
-- JOIN scan_file_tags t ON t.file_id = sf.id
-- WHERE t.tag IN ('jetTrkMaxFilter', 'WDecayFilter', 'pfCandCS')
-- GROUP BY sf.id
-- HAVING COUNT(DISTINCT t.tag) = 3;

-- 6. Most common tags across the whole catalog (what cuts show up most).
-- SELECT tag, COUNT(*) AS n_files
-- FROM scan_file_tags
-- GROUP BY tag
-- ORDER BY n_files DESC
-- LIMIT 15;
