use rusqlite::Connection;
use serde::{Deserialize, Serialize};
use std::fs;
use std::os::unix::fs::MetadataExt;
use std::sync::Mutex;
use tauri::Manager;
use std::path::PathBuf;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Paper {
    pub id: String,
    pub title: String,
    pub authors: String,
    pub year: i32,
    pub journal: String,
    pub abstract_text: String,
    pub doi: String,
    pub date_added: String,
    pub starred: i32,
    pub status: String,
    pub file_path: String,
    pub pages: i32,
    pub tags: String,
    pub category: String,
}

struct AppState {
    pdf_dir: PathBuf,
}

fn init_db(app: &tauri::App) -> Result<(), Box<dyn std::error::Error>> {
    let app_dir = app.path().app_data_dir().expect("无法获取应用数据目录");
    fs::create_dir_all(&app_dir)?;
    let db_path = app_dir.join("koan.db");
    let pdf_dir = app_dir.join("pdfs");
    fs::create_dir_all(&pdf_dir)?;

    let conn = Connection::open(&db_path)?;

    conn.execute_batch(
        "CREATE TABLE IF NOT EXISTS papers (
            id TEXT PRIMARY KEY,
            title TEXT NOT NULL,
            authors TEXT,
            year INTEGER,
            journal TEXT,
            abstract_text TEXT,
            doi TEXT,
            date_added TEXT,
            starred INTEGER DEFAULT 0,
            status TEXT DEFAULT 'unread',
            file_path TEXT DEFAULT '',
            pages INTEGER DEFAULT 0,
            tags TEXT DEFAULT '',
            category TEXT DEFAULT ''
        )"
    )?;

    app.manage(Mutex::new(conn));
    app.manage(AppState { pdf_dir });
    Ok(())
}

#[tauri::command]
fn get_papers(db: tauri::State<'_, Mutex<Connection>>) -> Result<Vec<Paper>, String> {
    let conn = db.lock().map_err(|e| e.to_string())?;
    let mut stmt = conn.prepare(
        "SELECT id, title, authors, year, journal, abstract_text, doi, date_added, 
         starred, status, file_path, pages, tags, category FROM papers ORDER BY date_added DESC"
    ).map_err(|e| e.to_string())?;

    let papers = stmt.query_map([], |row| {
        Ok(Paper {
            id: row.get(0)?,
            title: row.get(1)?,
            authors: row.get(2)?,
            year: row.get(3)?,
            journal: row.get(4)?,
            abstract_text: row.get(5)?,
            doi: row.get(6)?,
            date_added: row.get(7)?,
            starred: row.get(8)?,
            status: row.get(9)?,
            file_path: row.get(10)?,
            pages: row.get(11)?,
            tags: row.get(12)?,
            category: row.get(13)?,
        })
    }).map_err(|e| e.to_string())?.filter_map(|p| p.ok()).collect();

    Ok(papers)
}

#[tauri::command]
fn search_papers(db: tauri::State<'_, Mutex<Connection>>, keyword: String) -> Result<Vec<Paper>, String> {
    let conn = db.lock().map_err(|e| e.to_string())?;
    let search = format!("%{}%", keyword.to_lowercase());
    let mut stmt = conn.prepare(
        "SELECT id, title, authors, year, journal, abstract_text, doi, date_added, 
         starred, status, file_path, pages, tags, category FROM papers 
         WHERE LOWER(title) LIKE ?1 OR LOWER(authors) LIKE ?1 OR LOWER(journal) LIKE ?1 
         OR LOWER(tags) LIKE ?1 OR LOWER(category) LIKE ?1 OR LOWER(doi) LIKE ?1
         ORDER BY date_added DESC"
    ).map_err(|e| e.to_string())?;

    let papers = stmt.query_map([&search], |row| {
        Ok(Paper {
            id: row.get(0)?,
            title: row.get(1)?,
            authors: row.get(2)?,
            year: row.get(3)?,
            journal: row.get(4)?,
            abstract_text: row.get(5)?,
            doi: row.get(6)?,
            date_added: row.get(7)?,
            starred: row.get(8)?,
            status: row.get(9)?,
            file_path: row.get(10)?,
            pages: row.get(11)?,
            tags: row.get(12)?,
            category: row.get(13)?,
        })
    }).map_err(|e| e.to_string())?.filter_map(|p| p.ok()).collect();

    Ok(papers)
}

#[tauri::command]
fn add_paper_with_file(
    db: tauri::State<'_, Mutex<Connection>>,
    pdf_dir: tauri::State<'_, AppState>,
    paper: Paper,
    source_path: String,
) -> Result<(), String> {
    let mut final_paper = paper.clone();
    let mut pages = 0;
    let mut file_path = String::new();

    if !source_path.is_empty() {
        let dest_path = pdf_dir.pdf_dir.join(format!("{}.pdf", paper.id));
        if let Ok(_) = fs::copy(&source_path, &dest_path) {
            file_path = format!("{}.pdf", paper.id);
            if let Ok(metadata) = fs::metadata(&dest_path) {
                pages = std::cmp::max(1, (metadata.size() / 150000) as i32);
            }
        }
    }

    final_paper.file_path = file_path;
    final_paper.pages = pages;

    let conn = db.lock().map_err(|e| e.to_string())?;
    conn.execute(
        "INSERT INTO papers (id, title, authors, year, journal, abstract_text, doi, date_added, starred, status, file_path, pages, tags, category) 
         VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14)",
        rusqlite::params![
            final_paper.id, final_paper.title, final_paper.authors, final_paper.year, final_paper.journal,
            final_paper.abstract_text, final_paper.doi, final_paper.date_added, final_paper.starred,
            final_paper.status, final_paper.file_path, final_paper.pages, final_paper.tags, final_paper.category
        ],
    ).map_err(|e| e.to_string())?;
    Ok(())
}

#[tauri::command]
fn update_paper(db: tauri::State<'_, Mutex<Connection>>, paper: Paper) -> Result<(), String> {
    let conn = db.lock().map_err(|e| e.to_string())?;
    conn.execute(
        "UPDATE papers SET title=?2, authors=?3, year=?4, journal=?5, abstract_text=?6, doi=?7, 
         starred=?9, status=?10, file_path=?11, pages=?12, tags=?13, category=?14 WHERE id=?1",
        rusqlite::params![
            paper.id, paper.title, paper.authors, paper.year, paper.journal,
            paper.abstract_text, paper.doi, paper.date_added, paper.starred,
            paper.status, paper.file_path, paper.pages, paper.tags, paper.category
        ],
    ).map_err(|e| e.to_string())?;
    Ok(())
}

#[tauri::command]
fn delete_paper(db: tauri::State<'_, Mutex<Connection>>, id: String) -> Result<(), String> {
    let conn = db.lock().map_err(|e| e.to_string())?;
    conn.execute("DELETE FROM papers WHERE id=?1", [&id]).map_err(|e| e.to_string())?;
    Ok(())
}

#[tauri::command]
fn get_pdf_dir(pdf_dir: tauri::State<'_, AppState>) -> Result<String, String> {
    Ok(pdf_dir.pdf_dir.to_string_lossy().to_string())
}

#[tauri::command]
fn fetch_doi_info(doi: String) -> Result<Paper, String> {
    // 通过 CrossRef API 获取 DOI 信息
    let url = format!("https://api.crossref.org/works/{}", doi);
    
    let response = ureq::get(&url)
        .set("User-Agent", "KoanLiterature/1.0 (mailto:contact@koan.app)")
        .call()
        .map_err(|e| e.to_string())?;

    let body = response.into_string().map_err(|e| e.to_string())?;
    let json: serde_json::Value = serde_json::from_str(&body).map_err(|e| e.to_string())?;
    let message = &json["message"];
            
            let title = message["title"].as_array()
                .and_then(|arr| arr.first())
                .and_then(|v| v.as_str())
                .unwrap_or("")
                .to_string();
                
            let authors: Vec<String> = message["author"].as_array()
                .map(|arr| arr.iter()
                    .filter_map(|a| {
                        let given = a["given"].as_str().unwrap_or("");
                        let family = a["family"].as_str().unwrap_or("");
                        if family.is_empty() { None } else { Some(format!("{} {}", given, family).trim().to_string()) }
                    })
                    .collect())
                .unwrap_or_default();
                
            let year = message["published"].as_array()
                .and_then(|arr| arr.first())
                .and_then(|v| v["date-parts"].as_array())
                .and_then(|arr| arr.first())
                .and_then(|v| v.as_array())
                .and_then(|arr| arr.first())
                .and_then(|v| v.as_i64())
                .map(|y| y as i32)
                .unwrap_or(2024);
                
            let journal = message["container-title"].as_array()
                .and_then(|arr| arr.first())
                .and_then(|v| v.as_str())
                .unwrap_or("")
                .to_string();
                
            let abstract_text = message["abstract"].as_str()
                .unwrap_or("")
                .replace("<jats:p>", "")
                .replace("</jats:p>", "")
                .to_string();

            Ok(Paper {
                id: String::new(),
                title,
                authors: authors.join(", "),
                year,
                journal,
                abstract_text,
                doi: doi.clone(),
                date_added: String::new(),
                starred: 0,
                status: "unread".to_string(),
                file_path: String::new(),
                pages: 0,
                tags: String::new(),
                category: String::new(),
            })
}

fn main() {
    tauri::Builder::default()
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_dialog::init())
        .setup(|app| {
            init_db(app)?;
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![
            get_papers,
            search_papers,
            add_paper_with_file,
            update_paper,
            delete_paper,
            get_pdf_dir,
            fetch_doi_info
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}