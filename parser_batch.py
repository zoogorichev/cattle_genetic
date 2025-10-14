#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Пакетный парсер быков с обработкой в новых вкладках
"""

import csv
import time
import random
import os
import json
import re
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import (NoSuchElementException, 
                                       TimeoutException, 
                                       StaleElementReferenceException)

# Настройки Chrome
options = Options()
options.headless = False  # Для отладки
options.add_argument("--disable-blink-features=AutomationControlled")
options.add_argument("--window-size=1920,1080")
options.add_argument("--no-sandbox")
options.add_argument("--disable-dev-shm-usage")
options.add_argument("--disable-gpu")
options.add_argument("--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36")

# Конфигурационные параметры
MAX_PAGES = 1000  # Максимальное количество страниц
SAVE_INTERVAL = 10  # Сохранять прогресс каждые N страниц
PROFILE_DELAY = 0.5  # Пауза между обработкой профилей (секунды)

# Маркеры и требуемый порядок колонок
ORDERED_LOCI = [
    "TGLA227", "BM2113", "TGLA53", "ETH10", "SPS115", "TGLA122", "INRA23",
    "TGLA126", "BM1818", "ETH225", "BM1824", "CSRM60", "CSSM43", "ETH3",
    "ILST006", "HAUT27", "AMEL",
]

# Нормализация названий локусов
LOCUS_NORMALIZATION = {
    "INRA023": "INRA23",
    "ILSTS006": "ILST006",
    "SPS113": "SPS115",
}

# Регулярка для пары locus_allele1/allele2
PAIR_RE = re.compile(r"^([A-Za-z0-9]+)\s*[_\-]\s*([0-9]+)\s*/\s*([0-9]+)\s*$")

# Файлы для сохранения
csv_file = 'bulls_data.csv'
progress_file = 'progress.json'
links_file = 'bulls_links.json'  # Файл для сохранения всех ссылок

def normalize_locus(raw: str) -> str:
    """Нормализация названия локуса"""
    name = raw.strip().upper()
    name = LOCUS_NORMALIZATION.get(name, name)
    return name

def parse_profile_to_dict(profile_text: str) -> dict:
    """Разобрать строку профиля вида "BM1818_266/270, ..." -> {"BM1818": ("266","270"), ...}"""
    result = {}
    if not isinstance(profile_text, str) or not profile_text.strip():
        return result

    parts = [p.strip() for p in profile_text.split(",") if p.strip()]
    for part in parts:
        m = PAIR_RE.match(part)
        if not m:
            alt = re.sub(r"\s+", "_", part)
            m = PAIR_RE.match(alt)
        if not m:
            continue

        locus_raw, a1, a2 = m.group(1), m.group(2), m.group(3)
        locus = normalize_locus(locus_raw)
        if locus not in ORDERED_LOCI:
            continue
        result[locus] = (a1, a2)

    return result

def init_driver():
    """Инициализация драйвера с обработкой ошибок"""
    try:
        print("  Создание Chrome WebDriver...")
        driver = webdriver.Chrome(options=options)
        print("  WebDriver создан успешно")
        
        print("  Настройка таймаутов...")
        driver.set_page_load_timeout(60)
        driver.implicitly_wait(10)
        print("  Таймауты настроены")
        
        return driver
    except Exception as e:
        print(f"  Ошибка инициализации драйвера: {e}")
        print("  Возможные причины:")
        print("    - Chrome не установлен")
        print("    - ChromeDriver не найден в PATH")
        print("    - Недостаточно прав для запуска браузера")
        return None

def safe_get(driver, url, max_retries=3):
    """Безопасный переход на страницу"""
    for attempt in range(max_retries):
        try:
            print(f"  Переход на {url} (попытка {attempt + 1})")
            driver.get(url)
            print(f"  Страница загружена, проверяем заголовок...")
            
            # Проверяем, что страница загрузилась
            title = driver.title
            print(f"  Заголовок страницы: {title}")
            
            if "Быки России" in title or "Фильтры" in title:
                print(f"  ✓ Страница загружена успешно")
                return True
            else:
                print(f"  ⚠️ Неожиданный заголовок страницы")
                return True  # Все равно продолжаем
                
        except Exception as e:
            print(f"  Ошибка перехода (попытка {attempt + 1}): {e}")
            if attempt < max_retries - 1:
                print(f"  Ждем 3 секунды перед повторной попыткой...")
                time.sleep(3)
            else:
                print(f"  Все попытки исчерпаны")
                return False
    return False

def load_progress():
    """Загрузка прогресса"""
    if os.path.exists(progress_file):
        with open(progress_file, 'r', encoding='utf-8') as f:
            return json.load(f)
    return {'last_page': 0, 'processed_pages': [], 'collected_links': 0}

def save_progress(progress):
    """Сохранение прогресса"""
    with open(progress_file, 'w', encoding='utf-8') as f:
        json.dump(progress, f, ensure_ascii=False, indent=2)

def load_links():
    """Загрузка собранных ссылок"""
    if os.path.exists(links_file):
        with open(links_file, 'r', encoding='utf-8') as f:
            return json.load(f)
    return []

def save_links(links):
    """Сохранение собранных ссылок"""
    with open(links_file, 'w', encoding='utf-8') as f:
        json.dump(links, f, ensure_ascii=False, indent=2)

def append_to_csv(data_list):
    """Запись данных в CSV"""
    file_exists = os.path.isfile(csv_file)
    
    meta_keys = ['Идентификационный номер', 'Дата рождения', 'Ссылка']
    loci_keys = []
    for locus in ORDERED_LOCI:
        loci_keys.append(f"1_{locus}")
        loci_keys.append(f"2_{locus}")
    
    all_keys = meta_keys + loci_keys
    
    with open(csv_file, 'a', newline='', encoding='utf-8-sig') as output_file:
        writer = csv.DictWriter(output_file, all_keys, delimiter=';')
        if not file_exists:
            writer.writeheader()
        writer.writerows(data_list)

def collect_links_from_page(driver, page_num):
    """Сбор ссылок с одной страницы"""
    links = []
    try:
        print(f"  === СБОР ССЫЛОК СО СТРАНИЦЫ {page_num} ===")
        
        # Ждем загрузки AngularJS данных
        print("  Ожидаем загрузки данных...")
        WebDriverWait(driver, 30).until(
            EC.presence_of_element_located((By.XPATH, '//div[@ng-repeat="animal in animals"]'))
        )
        
        # Дополнительная проверка загрузки
        time.sleep(3)
        animal_rows = driver.find_elements(By.XPATH, '//div[@ng-repeat="animal in animals"]')
        print(f"  Найдено {len(animal_rows)} строк с быками")
        
        for i, row in enumerate(animal_rows):
            try:
                # Ищем ссылку на профиль
                inv_link = row.find_element(By.XPATH, './/a[contains(@ng-href, "/bulls/bull/")]')
                profile_url = inv_link.get_attribute('href')
                inv_number = inv_link.text.strip()
                
                # Ищем ID номер - упрощенный поиск
                id_number = 'Не найдено'
                try:
                    # Ищем по тексту в div элементах - ID обычно в формате US0018553781
                    all_divs = row.find_elements(By.XPATH, './/div')
                    for div in all_divs:
                        text = div.text.strip()
                        # Ищем ID по паттерну: 2-3 буквы + много цифр
                        if text and len(text) > 8 and any(c.isalpha() for c in text[:3]) and any(c.isdigit() for c in text[3:]):
                            id_number = text
                            break
                except:
                    pass
                
                # Ищем дату рождения - упрощенный поиск
                birth_date = 'Не найдено'
                try:
                    # Ищем дату в формате DD.MM.YYYY
                    all_divs = row.find_elements(By.XPATH, './/div')
                    for div in all_divs:
                        text = div.text.strip()
                        # Проверяем формат даты: DD.MM.YYYY
                        if text and '.' in text and len(text) == 10 and text.count('.') == 2:
                            parts = text.split('.')
                            if len(parts) == 3 and all(part.isdigit() for part in parts):
                                birth_date = text
                                break
                except:
                    pass
                
                links.append({
                    'url': profile_url,
                    'inv_number': inv_number,
                    'id_number': id_number,
                    'birth_date': birth_date,
                    'page': page_num
                })
                
                print(f"    {i+1}. Инв: {inv_number} | ID: {id_number} | Дата: {birth_date}")
                
                # Краткая диагностика только для первых 2 элементов
                if i < 2:
                    print(f"      Текст: {row.text[:100]}...")
                
            except Exception as e:
                print(f"    Ошибка обработки строки {i+1}: {e}")
                continue
        
        print(f"  === СОБРАНО {len(links)} ССЫЛОК СО СТРАНИЦЫ {page_num} ===")
        
    except Exception as e:
        print(f"Ошибка сбора ссылок со страницы {page_num}: {e}")
    
    return links

def process_profile_in_new_tab(driver, profile_info):
    """Обработка профиля в новой вкладке"""
    try:
        print(f"  Обрабатываем профиль: {profile_info['inv_number']}")
        
        # Открываем новую вкладку
        driver.execute_script("window.open('');")
        driver.switch_to.window(driver.window_handles[-1])
        
        # Переходим на профиль
        if not safe_get(driver, profile_info['url']):
            print(f"    ✗ Не удалось загрузить профиль")
            driver.close()
            driver.switch_to.window(driver.window_handles[0])
            return None
        
        time.sleep(2)  # Пауза для загрузки
        
        # Ищем микросателлитный профиль
        micro_profile = 'Не найдено'
        try:
            profile_row = driver.find_element(By.XPATH, '//td[contains(text(), "Микросателлитный профиль")]/following-sibling::td')
            micro_profile_text = profile_row.text.strip()
            
            if micro_profile_text and micro_profile_text != '':
                micro_profile = micro_profile_text
                print(f"    ✓ Найден микросателлитный профиль: {micro_profile[:100]}...")
            else:
                print(f"    ✗ Микросателлитный профиль пустой")
                driver.close()
                driver.switch_to.window(driver.window_handles[0])
                return None
                
        except Exception as e:
            print(f"    ✗ Микросателлитный профиль не найден: {e}")
            driver.close()
            driver.switch_to.window(driver.window_handles[0])
            return None
        
        # Парсим профиль
        parsed_profile = parse_profile_to_dict(micro_profile)
        
        if not parsed_profile:
            print(f"    ✗ Не удалось распарсить микросателлитный профиль")
            driver.close()
            driver.switch_to.window(driver.window_handles[0])
            return None
        
        print(f"    ✓ Парсинг профиля: найдено {len(parsed_profile)} локусов")
        
        # Создаем запись
        record = {
            'Идентификационный номер': profile_info['id_number'],
            'Дата рождения': profile_info['birth_date'],
            'Ссылка': profile_info['url']
        }
        
        # Добавляем колонки для каждого локуса
        for locus in ORDERED_LOCI:
            if locus in parsed_profile:
                record[f"1_{locus}"] = parsed_profile[locus][0]
                record[f"2_{locus}"] = parsed_profile[locus][1]
            else:
                record[f"1_{locus}"] = ""
                record[f"2_{locus}"] = ""
        
        # Закрываем вкладку и возвращаемся к основной
        driver.close()
        driver.switch_to.window(driver.window_handles[0])
        
        return record
        
    except Exception as e:
        print(f"Ошибка обработки профиля {profile_info['inv_number']}: {e}")
        try:
            driver.close()
            driver.switch_to.window(driver.window_handles[0])
        except:
            pass
        return None

def main():
    """Основная функция"""
    print("=== ПАКЕТНЫЙ ПАРСЕР БЫКОВ ===")
    
    # Инициализация
    print("Инициализация браузера...")
    driver = init_driver()
    if not driver:
        print("Не удалось инициализировать драйвер!")
        return
    
    base_url = 'https://xn--90aof1e.xn--p1ai/bulls/list'
    
    # Загружаем главную страницу
    print("Загрузка главной страницы...")
    if not safe_get(driver, base_url):
        print("Не удалось загрузить главную страницу!")
        driver.quit()
        return
    
    print("Главная страница загружена успешно!")
    
    # Загрузка прогресса
    progress = load_progress()
    all_links = load_links()
    
    print(f"Загружено {len(all_links)} ссылок из предыдущих сессий")
    
    # ЭТАП 1: Сбор всех ссылок
    print("\n=== ЭТАП 1: СБОР ВСЕХ ССЫЛОК ===")
    
    current_page = progress.get('last_page', 0) + 1
    processed_pages = set(progress.get('processed_pages', []))
    
    while current_page <= MAX_PAGES:
        print(f"\n--- Страница {current_page} ---")
        
        if current_page in processed_pages:
            print(f"  Страница {current_page} уже обработана, пропускаем")
            current_page += 1
            continue
        
        # Переходим на страницу
        if current_page == 1:
            page_url = base_url
        else:
            # Для AngularJS приложения используем JavaScript навигацию
            try:
                print(f"  Переход на страницу {current_page} через JavaScript...")
                driver.execute_script(f"goToPage({current_page})")
                time.sleep(3)  # Ждем загрузки
            except Exception as e:
                print(f"  Ошибка JavaScript навигации: {e}")
                
                # Fallback 1: попробуем кликнуть по кнопке "следующая страница"
                try:
                    print(f"  Пробуем кликнуть по кнопке 'следующая страница'...")
                    next_btn = driver.find_element(By.XPATH, '//a[@title="следующая страница" and contains(@ng-click, "goToPage")]')
                    driver.execute_script("arguments[0].click();", next_btn)
                    time.sleep(3)
                except Exception as e2:
                    print(f"  Ошибка клика по кнопке: {e2}")
                    
                    # Fallback 2: попробуем URL параметр
                    page_url = f"{base_url}?page={current_page}"
                    if not safe_get(driver, page_url):
                        print(f"  Не удалось загрузить страницу {current_page}")
                        break
                    continue
        
        # Проверяем, что мы на правильной странице
        if current_page > 1:
            try:
                # Проверяем, что данные изменились (не те же самые, что на предыдущей странице)
                first_row = driver.find_element(By.XPATH, '//div[@ng-repeat="animal in animals"][1]')
                first_text = first_row.text[:100]  # Берем первые 100 символов
                print(f"  Первая строка на странице {current_page}: {first_text}...")
            except Exception as e:
                print(f"  Ошибка проверки страницы: {e}")
        
        # Собираем ссылки
        page_links = collect_links_from_page(driver, current_page)
        
        if page_links:
            all_links.extend(page_links)
            save_links(all_links)
            
            # Обновляем прогресс
            processed_pages.add(current_page)
            progress['last_page'] = current_page
            progress['processed_pages'] = list(processed_pages)
            progress['collected_links'] = len(all_links)
            save_progress(progress)
            
            print(f"  Всего собрано ссылок: {len(all_links)}")
            current_page += 1
        else:
            print(f"  Не найдено ссылок на странице {current_page}, завершаем сбор")
            break
    
    print(f"\n=== СБОР ССЫЛОК ЗАВЕРШЕН ===")
    print(f"Всего собрано: {len(all_links)} ссылок")
    
    # ЭТАП 2: Обработка профилей
    print(f"\n=== ЭТАП 2: ОБРАБОТКА ПРОФИЛЕЙ ===")
    
    processed_count = 0
    successful_count = 0
    
    for i, profile_info in enumerate(all_links):
        # Показываем прогресс каждые 10 профилей или для первых 5
        if i < 5 or (i + 1) % 10 == 0:
            print(f"\nОбрабатываем профиль {i+1}/{len(all_links)} ({((i+1)/len(all_links)*100):.1f}%)")
        else:
            print(f"  {i+1}/{len(all_links)}", end=" ", flush=True)
        
        # Обрабатываем профиль в новой вкладке
        profile_data = process_profile_in_new_tab(driver, profile_info)
        
        if profile_data:
            # Сохраняем данные
            append_to_csv([profile_data])
            successful_count += 1
            if i < 5 or (i + 1) % 10 == 0:
                print(f"  ✓ Профиль сохранен")
        else:
            if i < 5 or (i + 1) % 10 == 0:
                print(f"  ✗ Профиль пропущен")
        
        processed_count += 1
        
        # Периодически сохраняем прогресс
        if processed_count % SAVE_INTERVAL == 0:
            print(f"\n  Обработано {processed_count}/{len(all_links)} профилей")
        
        # Пауза между профилями
        time.sleep(PROFILE_DELAY)
    
    # Завершение
    driver.quit()
    print(f"\n=== РАБОТА ЗАВЕРШЕНА ===")
    print(f"Обработано профилей: {processed_count}")
    print(f"Успешно сохранено: {successful_count}")
    print(f"Результаты сохранены в {csv_file}")

if __name__ == "__main__":
    main()
