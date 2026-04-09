#!/usr/bin/env python3
"""
Generate an interactive HTML index page for IGV reports.

This script scans the IGV reports directory and generates a sortable HTML page
that links to each report. It can optionally read a TSV file to mark flagged mutations.

TSV Format (optional):
    mutation_id    flagged    sample_id
    chr1-11243066-T-A    yes    SAMPLE_001
    chr1-150575252-G-A    no    SAMPLE_002
"""

import os
import csv
from pathlib import Path
from typing import Dict, List, Tuple


def parse_filename(filename: str) -> Dict[str, str]:
    """
    Parse an IGV report HTML filename into its components.
    
    Format: chr{number}-{position}-{ref}-{alt}-{annotation}-{gene}.html
    """
    # Remove .html extension
    name = filename[:-5]
    
    # Split on dashes, but we need to be careful because some fields might contain dashes
    parts = name.split('-')
    
    # Expected: chr#, position, ref, alt, annotation, gene
    if len(parts) >= 6:
        return {
            'chromosome': parts[0],
            'position': parts[1],
            'ref': parts[2],
            'alt': parts[3],
            'annotation': parts[4],
            'gene': parts[5],
            'filename': filename,
            'mutation_id': f"{parts[0]}-{parts[1]}-{parts[2]}-{parts[3]}"
        }
    return None


def split_sample_ids(sample_id_field: str) -> List[str]:
    """
    Split a sample_id field into one or more sample IDs.
    Accept both comma and semicolon separators.
    """
    if not sample_id_field:
        return []
    normalized = sample_id_field.replace(';', ',')
    return [part.strip() for part in normalized.split(',') if part.strip()]


def scan_maf_files(maf_dir: str) -> List[Dict[str, str]]:
    """
    Scan MAF files to extract individual mutation-sample_id combinations.
    Returns a list of dicts, each with 'mutation_id' and 'sample_id'
    """
    mutations = []
    
    if not os.path.isdir(maf_dir):
        print(f"Warning: MAF directory not found: {maf_dir}")
        return mutations
    
    for filename in sorted(os.listdir(maf_dir)):
        if not filename.endswith('.tsv'):
            continue
        
        # Extract mutation ID from filename (remove .tsv extension)
        # Filename format: chr1-11243066-T-A-intronic-MTOR.tsv
        # We want: chr1-11243066-T-A (stop at first dash after ref-alt)
        name_without_ext = filename[:-4]
        
        # Parse the mutation ID: chr-pos-ref-alt
        parts = name_without_ext.split('-')
        if len(parts) >= 4:
            mutation_id = f"{parts[0]}-{parts[1]}-{parts[2]}-{parts[3]}"
        else:
            mutation_id = name_without_ext  # fallback
        
        maf_path = os.path.join(maf_dir, filename)
        
        try:
            with open(maf_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    sample_id = row.get('sample_id', '').strip()
                    if sample_id:
                        mutations.append({
                            'mutation_id': mutation_id,
                            'sample_id': sample_id
                        })
        except Exception as e:
            print(f"Warning: Could not read MAF file {filename}: {e}")
            continue
    
    return mutations


def write_flagged_mutations(tsv_path: str, entries: List[Dict[str, str]]):
    """
    Write all mutation-sample_id entries back to the flagged_mutations.tsv file.
    Each entry represents one mutation-sample_id combination.
    """
    try:
        with open(tsv_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            
            # Write header
            writer.writerow(['mutation_id', 'flagged', 'sample_id'])
            
            # Write all entries, sorted by mutation_id then sample_id
            sorted_entries = sorted(entries, key=lambda x: (x['mutation_id'], x['sample_id']))
            
            for entry in sorted_entries:
                writer.writerow([
                    entry['mutation_id'],
                    entry['flagged'],
                    entry['sample_id']
                ])
        
        print(f"✅ Updated {tsv_path} with {len(entries)} mutation-sample_id combinations")
    except Exception as e:
        print(f"Error writing flagged mutations file: {e}")


def load_flagged_mutations(tsv_path: str, maf_entries: List[Dict[str, str]] = None) -> List[Dict[str, str]]:
    """
    Load flagged mutations from a TSV file.
    If maf_entries is provided, auto-populate with individual mutation-sample_id combinations from MAF files.
    Then merge with existing manual reviews from the TSV file.
    Returns a list of dicts, each with 'mutation_id', 'sample_id', and 'flagged' keys
    """
    entries = []
    entry_map = {}  # Track entries by (mutation_id, sample_id) for quick lookup
    
    # First, populate with all entries from MAF files
    if maf_entries:
        for entry in maf_entries:
            key = (entry['mutation_id'], entry['sample_id'])
            entries.append({
                'mutation_id': entry['mutation_id'],
                'sample_id': entry['sample_id'],
                'flagged': ''  # Empty status for auto-populated entries
            })
            entry_map[key] = entries[-1]  # Reference to the entry for updating
    
    if not os.path.exists(tsv_path):
        return entries
    
    try:
        with open(tsv_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                mutation_id = row.get('mutation_id', '').strip()
                flagged_status = row.get('flagged', '').strip().lower()
                sample_id = row.get('sample_id', '').strip()
                
                if not mutation_id or not sample_id:
                    continue
                
                # Normalize mutation ID to short format (chr-pos-ref-alt)
                parts = mutation_id.split('-')
                if len(parts) >= 4:
                    normalized_id = f"{parts[0]}-{parts[1]}-{parts[2]}-{parts[3]}"
                else:
                    normalized_id = mutation_id  # fallback
                
                key = (normalized_id, sample_id)
                
                # If this combination exists in our MAF data, update its flagged status
                if key in entry_map and flagged_status:
                    entry_map[key]['flagged'] = flagged_status
                # If it doesn't exist in MAF data but has a flagged status, add it
                elif flagged_status:
                    entries.append({
                        'mutation_id': normalized_id,
                        'sample_id': sample_id,
                        'flagged': flagged_status
                    })
                    entry_map[key] = entries[-1]
    except Exception as e:
        print(f"Warning: Could not read flagged mutations file: {e}")
    
    return entries


def get_igv_reports(reports_dir: str) -> List[Dict]:
    """
    Get all IGV report HTML files from the directory.
    """
    reports = []
    
    if not os.path.isdir(reports_dir):
        print(f"Error: Reports directory not found: {reports_dir}")
        return reports
    
    for filename in sorted(os.listdir(reports_dir)):
        if filename.endswith('.html'):
            parsed = parse_filename(filename)
            if parsed:
                reports.append(parsed)
    
    return reports


def generate_html(reports: List[Dict], flagged_mutations: Dict, output_path: str, base_url: str = ""):
    """
    Generate an interactive HTML page with a sortable table of IGV reports.
    """
    # Create table rows
    rows_html = []
    
    # Count samples per mutation_id
    sample_counts = {}
    for mutation_id, info in flagged_mutations.items():
        all_samples = info['flagged_samples'] | info['passed_samples'] | info['unreviewed_samples']
        sample_counts[mutation_id] = len(all_samples)
    
    for report in reports:
        mutation_id = report['mutation_id']
        info = flagged_mutations.get(mutation_id, {'status': '', 'sample_ids': []})
        flagged_status = info['status']
        flagged_samples = sorted(info['flagged_samples']) if info['flagged_samples'] else []
        passed_samples = sorted(info['passed_samples']) if info['passed_samples'] else []
        flagged_samples_str = ', '.join(flagged_samples) if flagged_samples else ''
        passed_samples_str = ', '.join(passed_samples) if passed_samples else ''
        sample_count = sample_counts.get(mutation_id, 0)
        
        # Determine row class and label based on review status
        if flagged_status == 'yes':
            row_class = 'flagged-row'
            flagged_label = '🚩 FLAGGED'
        elif flagged_status == 'no':
            row_class = 'passed-row'
            flagged_label = '✅ Passed'
        elif flagged_status == 'mixed':
            row_class = 'mixed-row'
            flagged_label = '⚠️ Mixed'
        else:
            row_class = ''
            flagged_label = 'pending'
        
        link = f"{base_url}{report['filename']}" if base_url else f"igv_reports/hotspots/{report['filename']}"
        
        rows_html.append(f"""
        <tr class="{row_class}">
            <td class="mutation-id">
                <div class="mutation-id-cell">
                    <a href="{link}" target="_blank">{mutation_id}</a>
                    <button class="copy-btn" type="button" data-mutation-id="{mutation_id}" aria-label="Copy mutation ID">
                        <svg viewBox="0 0 16 16" width="14" height="14" fill="currentColor" aria-hidden="true">
                            <path d="M5 1.5A1.5 1.5 0 0 1 6.5 0h5A1.5 1.5 0 0 1 13 1.5V2h.5A1.5 1.5 0 0 1 15 3.5v9A1.5 1.5 0 0 1 13.5 14H10v1.5A1.5 1.5 0 0 1 8.5 17h-5A1.5 1.5 0 0 1 2 15.5v-12A1.5 1.5 0 0 1 3.5 2H4v-.5z"/>
                            <path d="M4.5 3A.5.5 0 0 0 4 3.5V4h6v-.5a.5.5 0 0 0-.5-.5h-5z"/>
                        </svg>
                    </button>
                </div>
            </td>
            <td>{report['ref']} → {report['alt']}</td>
            <td>{report['annotation']}</td>
            <td>{report['gene']}</td>
            <td class="numeric">{sample_count}</td>
            <td class="flagged-status" data-sort="{flagged_status}">{flagged_label}</td>
            <td>{flagged_samples_str}</td>
            <td>{passed_samples_str}</td>
        </tr>
        """)
    
    rows_str = '\n'.join(rows_html)
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>IGV Reports Index</title>
    <style>
        * {{
            box-sizing: border-box;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
        }}
        
        h1 {{
            color: #333;
            margin-top: 0;
            border-bottom: 2px solid #007bff;
            padding-bottom: 10px;
        }}
        
        .info {{
            background-color: #f8f9fa;
            border-left: 4px solid #007bff;
            padding: 12px 16px;
            margin-bottom: 20px;
            border-radius: 4px;
            color: #666;
        }}
        
        .controls {{
            margin-bottom: 20px;
            display: flex;
            gap: 15px;
            align-items: center;
            flex-wrap: wrap;
        }}
        
        .search-box {{
            flex: 1;
            min-width: 250px;
        }}
        
        .search-box input {{
            width: 100%;
            padding: 10px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }}
        
        .search-box input:focus {{
            outline: none;
            border-color: #007bff;
            box-shadow: 0 0 0 3px rgba(0, 123, 255, 0.1);
        }}
        
        .filter-buttons {{
            display: flex;
            gap: 8px;
        }}
        
        button {{
            padding: 10px 16px;
            border: 1px solid #ddd;
            border-radius: 4px;
            background-color: white;
            cursor: pointer;
            font-size: 14px;
            transition: all 0.2s;
        }}
        
        button:hover {{
            border-color: #007bff;
            color: #007bff;
        }}
        
        button.active {{
            background-color: #007bff;
            color: white;
            border-color: #007bff;
        }}
        
        .stats {{
            font-size: 14px;
            color: #666;
            padding: 8px;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }}
        
        thead {{
            background-color: #f8f9fa;
            position: sticky;
            top: 0;
        }}
        
        th {{
            padding: 12px;
            text-align: left;
            border-bottom: 2px solid #dee2e6;
            font-weight: 600;
            color: #333;
            cursor: pointer;
            user-select: none;
            white-space: nowrap;
            position: relative;
        }}
        
        th:hover {{
            background-color: #e9ecef;
        }}
        
        th.sortable::after {{
            content: ' ↕';
            opacity: 0.3;
            margin-left: 4px;
        }}
        
        th.sort-asc::after {{
            content: ' ↑';
            opacity: 1;
            color: #007bff;
        }}
        
        th.sort-desc::after {{
            content: ' ↓';
            opacity: 1;
            color: #007bff;
        }}
        
        td {{
            padding: 11px 12px;
            border-bottom: 1px solid #dee2e6;
            font-size: 14px;
        }}
        
        tr:hover {{
            background-color: #f8f9fa;
        }}
        
        tr.flagged-row {{
            background-color: #fff3cd;
        }}
        
        tr.flagged-row:hover {{
            background-color: #ffe8a3;
        }}
        
        tr.passed-row {{
            background-color: #d4edda;
        }}
        
        tr.passed-row:hover {{
            background-color: #c3e6cb;
        }}
        
        tr.mixed-row {{
            background-color: #fff3cd;
            border-left: 4px solid #ffc107;
        }}
        
        tr.mixed-row:hover {{
            background-color: #ffe8a3;
        }}
        
        td.numeric {{
            text-align: right;
            font-variant-numeric: tabular-nums;
        }}
        
        a {{
            color: #007bff;
            text-decoration: none;
        }}
        
        a:hover {{
            text-decoration: underline;
        }}
        
        .mutation-id-cell {{
            display: flex;
            align-items: center;
            gap: 8px;
            flex-wrap: wrap;
        }}

        td.mutation-id a {{
            font-family: 'Monaco', 'Courier New', monospace;
            font-size: 13px;
        }}

        .copy-btn {{
            width: 26px;
            height: 26px;
            padding: 0;
            border: 1px solid #d8dee4;
            border-radius: 8px;
            background-color: #f6f8fa;
            color: #57606a;
            display: inline-flex;
            align-items: center;
            justify-content: center;
            cursor: pointer;
            transition: border-color 0.15s ease, background-color 0.15s ease, color 0.15s ease;
        }}

        .copy-btn:hover {{
            border-color: #959da5;
            background-color: #eef1f4;
            color: #24292f;
        }}

        .copy-btn.copied {{
            border-color: #28a745;
            background-color: #e6f4ea;
            color: #28a745;
        }}

        .sr-only {{
            position: absolute;
            width: 1px;
            height: 1px;
            padding: 0;
            margin: -1px;
            overflow: hidden;
            clip: rect(0, 0, 0, 0);
            white-space: nowrap;
            border: 0;
        }}

        td.flagged-status {{
            font-weight: 600;
        }}
        
        td.flagged-status:has-text('🚩') {{
            color: #d9534f;
        }}
        
        .flagged-row td.flagged-status {{
            color: #d9534f;
        }}
        
        .passed-row td.flagged-status {{
            color: #28a745;
        }}
        
        .hidden {{
            display: none;
        }}
        
        .footer {{
            margin-top: 30px;
            padding-top: 20px;
            border-top: 1px solid #dee2e6;
            font-size: 12px;
            color: #999;
        }}
        
        @media (max-width: 768px) {{
            .controls {{
                flex-direction: column;
            }}
            
            .search-box {{
                min-width: 100%;
            }}
            
            table {{
                font-size: 12px;
            }}
            
            th, td {{
                padding: 8px;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 IGV Reports Index</h1>
        
        <div class="info">
            Total variants: <strong id="total-count">{len(reports)}</strong> | 
            Flagged: <strong id="flagged-count">0</strong> | 
            Passed: <strong id="passed-count">0</strong>
        </div>
        
        <div class="controls">
            <div class="search-box">
                <input type="text" id="search-input" placeholder="Search by mutation ID, gene, chromosome...">
            </div>
            <div class="filter-buttons">
                <button id="filter-all" class="active">All</button>
                <button id="filter-flagged">Flagged Only</button>
            </div>
        </div>
        
        <div class="stats">
            Showing <span id="shown-count">{len(reports)}</span> of <span id="total-count-2">{len(reports)}</span>
        </div>
        
        <table id="reports-table">
            <thead>
                <tr>
                    <th class="sortable" data-column="mutation-id">Mutation ID</th>
                    <th class="sortable" data-column="change">Change</th>
                    <th class="sortable" data-column="annotation">Annotation</th>
                    <th class="sortable" data-column="gene">Gene</th>
                    <th class="sortable" data-column="sample-count">Samples</th>
                    <th class="sortable" data-column="flagged">Review Result</th>
                    <th class="sortable" data-column="flagged-samples">Flagged Samples</th>
                    <th class="sortable" data-column="passed-samples">Passed Samples</th>
                </tr>
            </thead>
            <tbody>
                {rows_str}
            </tbody>
        </table>
        
        <div class="footer">
            <p>Generated automatically by generate_igv_index.py</p>
            <p>Click on column headers to sort • Use search box to filter • Click filter buttons to show flagged variants</p>
        </div>
    </div>
    
    <script>
        const searchInput = document.getElementById('search-input');
        const filterAllBtn = document.getElementById('filter-all');
        const filterFlaggedBtn = document.getElementById('filter-flagged');
        const table = document.getElementById('reports-table');
        const tbody = table.querySelector('tbody');
        const rows = Array.from(tbody.querySelectorAll('tr'));
        
        let currentFilter = 'all';
        let sortColumn = null;
        let sortDirection = 'asc';
        
        // Count flagged mutations
        function updateStats() {{
            const visibleRows = rows.filter(r => !r.classList.contains('hidden'));
            const flaggedCount = rows.filter(r => r.classList.contains('flagged-row') || r.classList.contains('mixed-row')).length;
            const passedCount = rows.filter(r => r.classList.contains('passed-row')).length;
            document.getElementById('flagged-count').textContent = flaggedCount;
            document.getElementById('passed-count').textContent = passedCount;
            document.getElementById('shown-count').textContent = visibleRows.length;
        }}
        
        // Reset filter function
        function resetFilter() {{
            currentFilter = 'all';
            filterAllBtn.classList.add('active');
            filterFlaggedBtn.classList.remove('active');
            searchInput.value = '';
            applyFilters();
        }}
        
        // Filter and search functionality
        function applyFilters() {{
            const searchTerm = searchInput.value.toLowerCase();
            
            rows.forEach(row => {{
                let matches = true;
                
                // Search filter
                if (searchTerm) {{
                    const text = row.textContent.toLowerCase();
                    matches = text.includes(searchTerm);
                }}
                
                // Flagged filter
                if (currentFilter === 'flagged') {{
                    matches = matches && (row.classList.contains('flagged-row') || row.classList.contains('mixed-row'));
                }}
                
                row.classList.toggle('hidden', !matches);
            }});
            
            updateStats();
        }}
        
        // Sorting functionality
        function getSortValue(cell, column) {{
            const text = cell.textContent.trim();
            
            // Special handling for numeric columns
            if (column === 'position') {{
                return parseInt(text) || 0;
            }}
            
            // Review result column: sort flagged before passed
            if (column === 'flagged') {{
                if (text.includes('🚩')) return 0;
                if (text.includes('✅')) return 1;
                return 2;
            }}
            
            return text.toLowerCase();
        }}
        
        function sort(column) {{
            if (sortColumn === column) {{
                sortDirection = sortDirection === 'asc' ? 'desc' : 'asc';
            }} else {{
                sortColumn = column;
                sortDirection = 'asc';
            }}
            
            // Update header appearance
            document.querySelectorAll('th').forEach(th => {{
                th.classList.remove('sort-asc', 'sort-desc');
            }});
            
            const activeHeader = Array.from(document.querySelectorAll('th')).find(
                th => th.textContent.toLowerCase().includes(column.replace('-', ' ').split(': ')[1])
            );
            if (activeHeader) {{
                activeHeader.classList.add(`sort-${{sortDirection}}`);
            }}
            
            // Sort rows
            const visibleRows = rows.filter(r => !r.classList.contains('hidden'));
            visibleRows.sort((a, b) => {{
                let aValue, bValue;
                
                if (column === 'mutation-id') {{
                    aValue = a.cells[0].textContent.trim();
                    bValue = b.cells[0].textContent.trim();
                }} else if (column === 'change') {{
                    aValue = a.cells[1].textContent.trim();
                    bValue = b.cells[1].textContent.trim();
                }} else if (column === 'annotation') {{
                    aValue = a.cells[2].textContent.trim();
                    bValue = b.cells[2].textContent.trim();
                }} else if (column === 'gene') {{
                    aValue = a.cells[3].textContent.trim();
                    bValue = b.cells[3].textContent.trim();
                }} else if (column === 'sample-count') {{
                    aValue = parseInt(a.cells[4].textContent) || 0;
                    bValue = parseInt(b.cells[4].textContent) || 0;
                }} else if (column === 'flagged') {{
                    aValue = a.classList.contains('flagged-row') ? 0 : (a.classList.contains('mixed-row') ? 1 : (a.classList.contains('passed-row') ? 2 : 3));
                    bValue = b.classList.contains('flagged-row') ? 0 : (b.classList.contains('mixed-row') ? 1 : (b.classList.contains('passed-row') ? 2 : 3));
                }} else if (column === 'flagged-samples') {{
                    aValue = a.cells[6].textContent.trim();
                    bValue = b.cells[6].textContent.trim();
                }} else if (column === 'passed-samples') {{
                    aValue = a.cells[7].textContent.trim();
                    bValue = b.cells[7].textContent.trim();
                }}
                
                if (aValue < bValue) return sortDirection === 'asc' ? -1 : 1;
                if (aValue > bValue) return sortDirection === 'asc' ? 1 : -1;
                return 0;
            }});
            
            // Re-insert sorted rows
            visibleRows.forEach(row => tbody.appendChild(row));
        }}
        
        // Event listeners
        document.querySelectorAll('th.sortable').forEach(th => {{
            th.addEventListener('click', () => {{
                sort(th.dataset.column);
            }});
        }});
        
        searchInput.addEventListener('input', applyFilters);
        
        filterAllBtn.addEventListener('click', () => {{
            currentFilter = 'all';
            filterAllBtn.classList.add('active');
            filterFlaggedBtn.classList.remove('active');
            applyFilters();
        }});
        
        filterFlaggedBtn.addEventListener('click', () => {{
            currentFilter = 'flagged';
            filterFlaggedBtn.classList.add('active');
            filterAllBtn.classList.remove('active');
            applyFilters();
        }});

        function copyToClipboard(text) {{
            if (navigator.clipboard && navigator.clipboard.writeText) {{
                return navigator.clipboard.writeText(text);
            }}
            const textarea = document.createElement('textarea');
            textarea.value = text;
            textarea.style.position = 'fixed';
            textarea.style.top = '0';
            textarea.style.left = '0';
            textarea.style.opacity = '0';
            document.body.appendChild(textarea);
            textarea.focus();
            textarea.select();
            return new Promise((resolve, reject) => {{
                try {{
                    document.execCommand('copy');
                    resolve();
                }} catch (err) {{
                    reject(err);
                }} finally {{
                    document.body.removeChild(textarea);
                }}
            }});
        }}

        document.querySelectorAll('.copy-btn').forEach(button => {{
            button.addEventListener('click', () => {{
                const mutationId = button.dataset.mutationId;
                copyToClipboard(mutationId)
                    .then(() => {{
                        button.textContent = 'Copied';
                        button.classList.add('copied');
                        setTimeout(() => {{
                            button.textContent = 'Copy';
                            button.classList.remove('copied');
                        }}, 1200);
                    }})
                    .catch(() => {{
                        alert('Unable to copy Mutation ID to clipboard.');
                    }});
            }});
        }});
        
        // Initialize stats
        updateStats();
    </script>
</body>
</html>
"""
    
    with open(output_path, 'w') as f:
        f.write(html_content)
    
    flagged_count = sum(1 for info in flagged_mutations.values() if info.get('status') == 'yes')
    passed_count = sum(1 for info in flagged_mutations.values() if info.get('status') == 'no')

    print(f"✅ Generated index page: {output_path}")
    print(f"   - {len(reports)} variants indexed")
    print(f"   - {flagged_count} flagged mutations")
    print(f"   - {passed_count} passed mutations")


def entries_to_dict(entries: List[Dict[str, str]]) -> Dict[str, Dict[str, object]]:
    """
    Convert list of entries to dict format for HTML generation.
    Groups by mutation_id and collects sample_ids separated by flagged status.
    Determines overall status: mixed, yes, no, or empty.
    """
    mutations = {}
    
    for entry in entries:
        mutation_id = entry['mutation_id']
        sample_id = entry['sample_id']
        flagged = entry['flagged']
        
        if mutation_id not in mutations:
            mutations[mutation_id] = {
                'statuses': set(),
                'flagged_samples': set(),
                'passed_samples': set(),
                'unreviewed_samples': set()
            }
        
        # Track all statuses for this mutation
        if flagged:
            mutations[mutation_id]['statuses'].add(flagged)
        
        # Categorize samples by flagged status
        if flagged == 'yes':
            mutations[mutation_id]['flagged_samples'].add(sample_id)
        elif flagged == 'no':
            mutations[mutation_id]['passed_samples'].add(sample_id)
        else:
            mutations[mutation_id]['unreviewed_samples'].add(sample_id)
    
    # Determine overall status for each mutation
    for mutation_id, info in mutations.items():
        statuses = info['statuses']
        if len(statuses) > 1:
            info['status'] = 'mixed'
        elif len(statuses) == 1:
            info['status'] = next(iter(statuses))
        else:
            info['status'] = ''
        
        # Create combined sample_ids for backward compatibility
        all_samples = info['flagged_samples'] | info['passed_samples'] | info['unreviewed_samples']
        info['sample_ids'] = all_samples
    
    return mutations


def main():
    # Configuration
    REPORTS_DIR = "docs/igv_reports/hotspots"
    MAF_DIR = "data/grouped_maf"
    TSV_FILE = "flagged_mutations.tsv"
    OUTPUT_FILE = "docs/index.html"
    BASE_URL = ""  # For GitHub Pages, use: "https://github.com/username/panea/blob/main/docs/igv_reports/"
    
    # Load data
    print(f"📂 Scanning reports in {REPORTS_DIR}...")
    reports = get_igv_reports(REPORTS_DIR)
    print(f"   Found {len(reports)} IGV reports")
    
    print(f"📊 Scanning MAF files in {MAF_DIR}...")
    maf_entries = scan_maf_files(MAF_DIR)
    print(f"   Found {len(maf_entries)} mutation-sample_id combinations in MAF files")
    
    print(f"📄 Loading flagged mutations from {TSV_FILE}...")
    flagged_entries = load_flagged_mutations(TSV_FILE, maf_entries)
    if flagged_entries:
        print(f"   Loaded {len(flagged_entries)} total mutation-sample_id combinations")
    else:
        print(f"   No entries found")
    
    # Write back the complete list to flagged_mutations.tsv
    print(f"💾 Updating {TSV_FILE} with all mutation-sample_id combinations...")
    write_flagged_mutations(TSV_FILE, flagged_entries)
    
    # Generate HTML
    print(f"🔨 Generating HTML index...")
    flagged_mutations_dict = entries_to_dict(flagged_entries)
    generate_html(reports, flagged_mutations_dict, OUTPUT_FILE, BASE_URL)
    print("✅ Done!")


if __name__ == "__main__":
    main()
