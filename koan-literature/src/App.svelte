<script lang="ts">
 import { onMount } from "svelte";
 import { invoke } from "@tauri-apps/api/core";
 import { open } from "@tauri-apps/plugin-dialog";

 interface Paper {
  id: string;
  title: string;
  authors: string;
  year: number;
  journal: string;
  abstract_text: string;
  doi: string;
  date_added: string;
  starred: number;
  status: string;
  file_path: string;
  pages: number;
  tags: string;
  category: string;
 }

 let papers: Paper[] = $state([]);
 let formMode: 'add' | 'edit' = $state('add');
 let pdfDir = $state("");
 let searchKeyword = $state("");
 let showSearch = $state(false);
 
 // 表单字段
 let formTitle = $state("");
 let formAuthors = $state("");
 let formYear = $state(2024);
 let formJournal = $state("");
 let formAbstract = $state("");
 let formDoi = $state("");
 let formStatus = $state("unread");
 let formTags = $state("");
 let formCategory = $state("");
 let selectedFile = $state("");
 let editingId = $state("");

 // 分类选项
 const categories = ["", "眼科", "神经科", "精神科", "肿瘤", "心血管", "代谢", "其他"];

 async function loadPapers() {
  try {
   papers = await invoke<Paper[]>("get_papers");
   pdfDir = await invoke<string>("get_pdf_dir");
  } catch (error) {
   console.error("加载失败:", error);
  }
 }

 async function handleSearch() {
  if (!searchKeyword.trim()) {
   await loadPapers();
   return;
  }
  try {
   papers = await invoke<Paper[]>("search_papers", { keyword: searchKeyword });
  } catch (error) {
   console.error("搜索失败:", error);
  }
 }

 async function handleFileSelect() {
  try {
   const selected = await open({
    multiple: false,
    filters: [{ name: 'PDF', extensions: ['pdf'] }]
   });
   if (selected) {
    selectedFile = selected as string;
   }
  } catch (error) {
   console.error("选择文件失败:", error);
  }
 }

 async function fetchDoiInfo() {
  if (!formDoi.trim()) return;
  try {
   const paper = await invoke<Paper>("fetch_doi_info", { doi: formDoi });
   formTitle = paper.title || formTitle;
   formAuthors = paper.authors || formAuthors;
   formYear = paper.year || formYear;
   formJournal = paper.journal || formJournal;
   formAbstract = paper.abstract_text || formAbstract;
  } catch (error) {
   console.error("获取DOI信息失败:", error);
  }
 }

 function openAddForm() {
  formMode = 'add';
  formTitle = "";
  formAuthors = "";
  formYear = 2024;
  formJournal = "";
  formAbstract = "";
  formDoi = "";
  formStatus = "unread";
  formTags = "";
  formCategory = "";
  selectedFile = "";
  editingId = "";
 }

 function openEditForm(paper: Paper) {
  formMode = 'edit';
  formTitle = paper.title;
  formAuthors = paper.authors;
  formYear = paper.year;
  formJournal = paper.journal;
  formAbstract = paper.abstract_text;
  formDoi = paper.doi;
  formStatus = paper.status;
  formTags = paper.tags;
  formCategory = paper.category;
  selectedFile = "";
  editingId = paper.id;
 }

 async function handleSubmit() {
  if (!formTitle.trim()) return;

  const paper: Paper = {
   id: editingId || crypto.randomUUID(),
   title: formTitle,
   authors: formAuthors || "未知作者",
   year: formYear,
   journal: formJournal,
   abstract_text: formAbstract,
   doi: formDoi,
   date_added: new Date().toISOString(),
   starred: 0,
   status: formStatus,
   file_path: "",
   pages: 0,
   tags: formTags,
   category: formCategory
  };

  try {
   if (formMode === 'add') {
    await invoke("add_paper_with_file", { paper, sourcePath: selectedFile });
   } else {
    const existing = papers.find(p => p.id === editingId);
    paper.file_path = existing?.file_path || "";
    paper.pages = existing?.pages || 0;
    await invoke("update_paper", { paper });
   }
   openAddForm();
   await loadPapers();
  } catch (error) {
   console.error("保存失败:", error);
  }
 }

 async function handleDelete(id: string) {
  if (!confirm("确定删除这篇文献吗？")) return;
  try {
   await invoke("delete_paper", { id });
   await loadPapers();
  } catch (error) {
   console.error("删除失败:", error);
  }
 }

 async function toggleStar(paper: Paper) {
  paper.starred = paper.starred ? 0 : 1;
  try {
   await invoke("update_paper", { paper });
   await loadPapers();
  } catch (error) {
   console.error("更新失败:", error);
  }
 }

 function getStatusText(s: string): string {
  return s === 'unread' ? '未读' : s === 'reading' ? '在读' : '已读';
 }

 function getStatusColor(s: string): string {
  return s === 'unread' ? 'bg-slate-100 text-slate-600' : s === 'reading' ? 'bg-amber-100 text-amber-700' : 'bg-green-100 text-green-700';
 }

 onMount(() => {
  loadPapers();
 });
</script>

<main class="min-h-screen bg-slate-50 p-8 font-sans">
 <div class="max-w-3xl mx-auto">
  <div class="flex items-center justify-between mb-6">
   <h1 class="text-2xl font-bold text-slate-800 flex items-center gap-2">
    <span class="text-amber-500">⚓</span> 靠岸文献
   </h1>
   <div class="flex items-center gap-2">
    <button onclick={() => { showSearch = !showSearch; if(!showSearch) searchKeyword = ''; loadPapers(); }} 
     class="bg-slate-200 text-slate-700 px-3 py-2 rounded-md text-sm hover:bg-slate-300">
     {showSearch ? '🔍' : '🔍'} 搜索
    </button>
    <button onclick={openAddForm} class="bg-amber-600 text-white px-4 py-2 rounded-md text-sm font-medium hover:bg-amber-700">
     + 添加文献
    </button>
   </div>
  </div>

  <!-- 搜索栏 -->
  {#if showSearch}
   <div class="bg-white p-4 rounded-lg shadow-sm border border-slate-200 mb-6 flex gap-2">
    <input bind:value={searchKeyword} type="text" placeholder="搜索标题/作者/期刊/标签/分类/DOI..." 
     class="flex-1 border border-slate-300 rounded-md px-3 py-2 text-sm"
     onkeydown={(e) => e.key === 'Enter' && handleSearch()} />
    <button onclick={handleSearch} class="bg-blue-600 text-white px-4 py-2 rounded-md text-sm hover:bg-blue-700">
     搜索
    </button>
   </div>
  {/if}

  <!-- 表单 -->
  <div class="bg-white p-6 rounded-lg shadow-sm border border-slate-200 mb-6">
   <h2 class="text-lg font-semibold mb-4">{formMode === 'add' ? '添加新文献' : '编辑文献'}</h2>
   <div class="grid grid-cols-2 gap-4">
    <div class="col-span-2">
     <label class="block text-sm text-slate-600 mb-1">标题 *</label>
     <input bind:value={formTitle} type="text" class="w-full border border-slate-300 rounded-md px-3 py-2 text-sm" />
    </div>
    <div>
     <label class="block text-sm text-slate-600 mb-1">作者</label>
     <input bind:value={formAuthors} type="text" class="w-full border border-slate-300 rounded-md px-3 py-2 text-sm" />
    </div>
    <div class="flex gap-2">
     <div class="flex-1">
      <label class="block text-sm text-slate-600 mb-1">年份</label>
      <input bind:value={formYear} type="number" class="w-full border border-slate-300 rounded-md px-3 py-2 text-sm" />
     </div>
     <div class="flex-1">
      <label class="block text-sm text-slate-600 mb-1">DOI</label>
      <input bind:value={formDoi} type="text" placeholder="10.xxxx/xxxxx" class="w-full border border-slate-300 rounded-md px-3 py-2 text-sm" />
     </div>
    </div>
    <button onclick={fetchDoiInfo} class="text-sm text-blue-600 hover:text-blue-700 h-8 self-end -mt-2">
     📥 从DOI自动填充
    </button>
    <div>
     <label class="block text-sm text-slate-600 mb-1">期刊</label>
     <input bind:value={formJournal} type="text" class="w-full border border-slate-300 rounded-md px-3 py-2 text-sm" />
    </div>
    <div>
     <label class="block text-sm text-slate-600 mb-1">分类</label>
     <select bind:value={formCategory} class="w-full border border-slate-300 rounded-md px-3 py-2 text-sm">
      {#each categories as cat}
       <option value={cat}>{cat || '选择分类'}</option>
      {/each}
     </select>
    </div>
    <div>
     <label class="block text-sm text-slate-600 mb-1">状态</label>
     <select bind:value={formStatus} class="w-full border border-slate-300 rounded-md px-3 py-2 text-sm">
      <option value="unread">未读</option>
      <option value="reading">在读</option>
      <option value="read">已读</option>
     </select>
    </div>
    <div>
     <label class="block text-sm text-slate-600 mb-1">标签 (逗号分隔)</label>
     <input bind:value={formTags} type="text" placeholder="如: OCT, 抑郁症, 视网膜" class="w-full border border-slate-300 rounded-md px-3 py-2 text-sm" />
    </div>
    <div>
     <label class="block text-sm text-slate-600 mb-1">PDF文件</label>
     <button onclick={handleFileSelect} class="w-full border border-slate-300 rounded-md px-3 py-2 text-sm text-left hover:bg-slate-50 truncate">
      {selectedFile ? selectedFile.split(/[/\\]/).pop() : '点击选择PDF文件'}
     </button>
    </div>
    <div class="col-span-2">
     <label class="block text-sm text-slate-600 mb-1">摘要</label>
     <textarea bind:value={formAbstract} rows="3" class="w-full border border-slate-300 rounded-md px-3 py-2 text-sm"></textarea>
    </div>
   </div>
   <div class="flex gap-3 mt-4">
    <button onclick={handleSubmit} class="bg-amber-600 text-white px-4 py-2 rounded-md text-sm font-medium hover:bg-amber-700">
     {formMode === 'add' ? '添加文献' : '保存修改'}
    </button>
    <button onclick={openAddForm} class="bg-slate-200 text-slate-700 px-4 py-2 rounded-md text-sm font-medium hover:bg-slate-300">
     取消
    </button>
   </div>
  </div>

  <!-- 文献列表 -->
  <div class="space-y-3">
   {#each papers as paper}
    <div class="bg-white p-4 rounded-lg shadow-sm border border-slate-200 hover:border-amber-300 transition-colors">
     <div class="flex items-start justify-between">
      <div class="flex-1">
       <h3 class="font-semibold text-slate-800">{paper.title}</h3>
       <p class="text-xs text-slate-500 mt-1">{paper.authors} · {paper.year} · {paper.journal}</p>
       <div class="flex items-center gap-2 mt-2 flex-wrap">
        <span class="text-xs px-2 py-0.5 rounded {getStatusColor(paper.status)}">{getStatusText(paper.status)}</span>
        {#if paper.file_path}<span class="text-xs bg-blue-100 text-blue-700 px-2 py-0.5 rounded">📄 PDF</span>{/if}
        {#if paper.doi}<span class="text-xs bg-purple-100 text-purple-700 px-2 py-0.5 rounded">DOI</span>{/if}
        {#if paper.category}<span class="text-xs bg-cyan-100 text-cyan-700 px-2 py-0.5 rounded">{paper.category}</span>{/if}
        {#if paper.tags}
         {#each paper.tags.split(',') as tag}
          <span class="text-xs bg-green-100 text-green-700 px-2 py-0.5 rounded">{tag.trim()}</span>
         {/each}
        {/if}
       </div>
      </div>
      <div class="flex items-center gap-2 ml-4">
       <button onclick={() => toggleStar(paper)} class="text-lg">{paper.starred ? '⭐' : '☆'}</button>
       <button onclick={() => openEditForm(paper)} class="text-sm text-slate-500 hover:text-blue-600">编辑</button>
       <button onclick={() => handleDelete(paper.id)} class="text-sm text-slate-500 hover:text-red-600">删除</button>
      </div>
     </div>
    </div>
   {:else}
    <p class="text-center text-slate-400 text-sm py-10">暂无文献，请添加第一篇</p>
   {/each}
  </div>
 </div>
</main>